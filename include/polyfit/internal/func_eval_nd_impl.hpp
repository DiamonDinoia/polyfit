#pragma once

#include "impl_common.hpp"

// Out-of-class definitions for poly_eval::FuncEvalND.

namespace poly_eval {

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : coeffsFlat(), coeffsMd{coeffsFlat.data(), Extents{}} {
    static_assert(takesNdInput_v<Func>, "FuncEvalND requires fixed-size indexable ND input and output types");
    initialize(detail::CompileTimeCountTag{}, f, a, b);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b)
    : coeffsFlat(storageRequired(detail::validatePositiveCoeffCount(nCoeffsPerAxis))),
      coeffsMd{coeffsFlat.data(), makeExtents(nCoeffsPerAxis)} {
    static_assert(takesNdInput_v<Func>, "FuncEvalND requires fixed-size indexable ND input and output types");
    initialize(detail::RuntimeCountTag{}, f, nCoeffsPerAxis, a, b);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::initialize(detail::CompileTimeCountTag, Func f, const InputType &a,
                                                                  const InputType &b) {
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    detail::validateDomain(a, b);
    DomainParams dp;
    computeScaling(a, b, dp);
    buildCoeffs(static_cast<int>(NCOEFFS), f, dp);
    // Only the Always mode applies fusion up front; Auto keeps the baseline
    // eval-time mapFromDomain path so multi-panel fits don't pay the per-leaf
    // fusion cost by default.
    if constexpr (FUSION_MODE == FusionMode::Always) {
        fuseNDDomain(dp, static_cast<int>(NCOEFFS));
    }
    if constexpr (kStoresDomain) domain_ = dp;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::initialize(detail::RuntimeCountTag, Func f, int nCoeffsPerAxis,
                                                                  const InputType &a, const InputType &b) {
    static_assert(NCOEFFS == 0, "Runtime coefficient count is only valid for runtime-sized evaluators");
    detail::validateDomain(a, b);
    DomainParams dp;
    computeScaling(a, b, dp);
    buildCoeffs(nCoeffsPerAxis, f, dp);
    if constexpr (FUSION_MODE == FusionMode::Always) {
        fuseNDDomain(dp, nCoeffsPerAxis);
    }
    if constexpr (kStoresDomain) domain_ = dp;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr std::size_t FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::nCoeffsPerAxis() const noexcept {
    return static_cast<std::size_t>(coeffsMd.extent(0));
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::FuncEvalND(const FuncEvalND &other)
    : domain_(other.domain_), coeffsFlat(other.coeffsFlat),
      coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        domain_ = other.domain_;
        coeffsFlat = other.coeffsFlat;
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::FuncEvalND(FuncEvalND &&other) noexcept
    : domain_(std::move(other.domain_)),
      coeffsFlat(std::move(other.coeffsFlat)), coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        domain_ = std::move(other.domain_);
        coeffsFlat = std::move(other.coeffsFlat);
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(
    const InputType &x) const {
    return evalPoint<SIMD, InputType>(x);
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<bool SIMD, class Point, class>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(const Point &x) const {
    return evalPoint<SIMD>(x);
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<class... Coords, class>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(Coords... coords) const {
    return fromCanonicalOutput(evalCanonical(CanonicalInput{static_cast<InputScalar>(coords)...}));
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<bool SIMD, class Point>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::OutputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::evalPoint(const Point &x) const {
    return fromCanonicalOutput(evalCanonical<SIMD>(toCanonicalInput(x)));
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::CanonicalOutput
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::evalCanonical(const CanonicalInput &x) const noexcept {
    const int nCoeffsRt = (NCOEFFS ? static_cast<int>(NCOEFFS) : static_cast<int>(coeffsMd.extent(0)));
    const auto xm = mapFromDomain(x);
    if constexpr (SK == ScalarKernel::Horner) {
        return poly_eval::horner<NCOEFFS, SIMD, CanonicalOutput>(xm, coeffsMd, nCoeffsRt);
    } else {
        // Hybrid per axis. hybrid_nd internally falls back to horner for
        // runtime NCOEFFS or NCOEFFS ≤ 4 where Estrin has no critical-path
        // budget to spend.
        return poly_eval::hybrid_nd<NCOEFFS, CanonicalOutput, CanonicalInput, Mdspan, HK>(
            xm, coeffsMd, nCoeffsRt);
    }
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(const CanonicalInput *pts, CanonicalOutput *out,
                                                     std::size_t count) const noexcept {
    // Across-points vectorization: only profitable when OUT_DIM == 1, where
    // the existing OUT-dim SIMD path degenerates to scalar FMAs. For wider
    // outputs the per-point evaluator already vectorises across OUT_DIM.
    using batch_t = xsimd::batch<Scalar>;
    constexpr std::size_t B = batch_t::size;
    if constexpr (OUT_DIM != 1 || B <= 1) {
        for (std::size_t i = 0; i < count; ++i) out[i] = evalCanonical<>(pts[i]);
        return;
    } else {
        constexpr std::size_t kAlign = batch_t::arch_type::alignment();
        // Unroll factor across batches: independent Horner accumulator
        // chains per outer iteration so the FMA-throughput-bound inner
        // loop can run multiple chains in parallel and increase ILP.
        // Higher Dim builds deeper nested accumulator state per chain,
        // so unrolling spills registers; restrict unroll to low Dim.
        constexpr std::size_t U = (DIM <= 2) ? 2 : 1;
        const int nCoeffsRt = nCoeffsPerAxis();

        auto loadPointsAt = [&](std::size_t base) PF_ALWAYS_INLINE_LAMBDA {
            alignas(kAlign) Scalar buf[DIM][B];
            for (std::size_t b = 0; b < B; ++b)
                for (std::size_t d = 0; d < DIM; ++d)
                    buf[d][b] = pts[base + b][d];

            std::array<batch_t, DIM> x_v;
            for (std::size_t d = 0; d < DIM; ++d) {
                batch_t v = batch_t::load_aligned(buf[d]);
                if constexpr (FUSION_MODE != FusionMode::Always) {
                    if (!domain_.identityDomain) {
                        v = polyfit::internal::helpers::mapFromDomainScalar<batch_t, Scalar>(
                            v, domain_.invSpan[d], domain_.sumEndpoints[d]);
                    }
                }
                x_v[d] = v;
            }
            return x_v;
        };

        auto storeBatch = [&](std::size_t base, batch_t res) PF_ALWAYS_INLINE_LAMBDA {
            alignas(kAlign) Scalar outbuf[B];
            res.store_aligned(outbuf);
            for (std::size_t b = 0; b < B; ++b) out[base + b][0] = outbuf[b];
        };

        std::size_t i = 0;
        for (; i + U * B <= count; i += U * B) {
            std::array<std::array<batch_t, DIM>, U> x_vU;
            for (std::size_t u = 0; u < U; ++u) x_vU[u] = loadPointsAt(i + u * B);

            std::array<batch_t, U> resU;
            for (std::size_t u = 0; u < U; ++u) {
                resU[u] = detail::horner_nd_acrossPts<DIM, NCOEFFS, batch_t>(
                    x_vU[u], coeffsMd, nCoeffsRt);
            }

            for (std::size_t u = 0; u < U; ++u) storeBatch(i + u * B, resU[u]);
        }
        for (; i + B <= count; i += B) {
            auto x_v = loadPointsAt(i);
            batch_t res = detail::horner_nd_acrossPts<DIM, NCOEFFS, batch_t>(
                x_v, coeffsMd, nCoeffsRt);
            storeBatch(i, res);
        }
        // Tail: scalar fallback for the leftover < B points.
        for (; i < count; ++i) out[i] = evalCanonical<>(pts[i]);
    }
}

#if defined(__cpp_lib_span) && (__cpp_lib_span >= 202002L)
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(std::span<const CanonicalInput> pts,
                                                     std::span<CanonicalOutput> out) const {
    if (pts.size() != out.size()) {
        throw std::invalid_argument("Input and output spans must have equal length");
    }
    operator()(pts.data(), out.data(), pts.size());
}
#endif

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<class Points, class Outputs, class>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::operator()(const Points &pts, Outputs &out) const {
    if (pts.size() != out.size()) {
        throw std::invalid_argument("Input and output containers must have equal length");
    }

    const auto *ptsData = detail::FixedContainerAccess<Points>::data(pts);
    auto *outData = detail::FixedContainerAccess<Outputs>::data(out);
    for (std::size_t i = 0; i < static_cast<std::size_t>(pts.size()); ++i) {
        outData[i] = detail::fixedContainerCast<detail::data_value_t<Outputs>>(
            evalCanonical(toCanonicalInput(ptsData[i])));
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::coeffImpl(
    const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffsMd[std::array<std::size_t, DIM + 1>{static_cast<std::size_t>(idx[I])..., k}];
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::coeff(
    const IdxArray &idx, std::size_t k) noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<DIM>{});
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<typename IdxArray, std::size_t... I>
constexpr const typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::coeffImpl(
    const IdxArray &idx, std::size_t k, std::index_sequence<I...>) const noexcept {
    return coeffsMd[std::array<std::size_t, DIM + 1>{static_cast<std::size_t>(idx[I])..., k}];
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<class IdxArray>
[[nodiscard]] constexpr const typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::coeff(
    const IdxArray &idx, std::size_t k) const noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<DIM>{});
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::makeExtents(int nCoeffsPerAxis) noexcept -> Extents {
    if constexpr (IS_STATIC) {
        return detail::makeStaticExtents<NCOEFFS, DIM, OUT_DIM>(std::make_index_sequence<DIM>{});
    } else {
        return makeExtents(nCoeffsPerAxis, std::make_index_sequence<DIM + 1>{});
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<std::size_t... Is>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept -> Extents {
    return Extents{(Is < DIM ? static_cast<std::size_t>(nCoeffsPerAxis) : static_cast<std::size_t>(OUT_DIM))...};
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr std::size_t FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::storageRequired(const int nCoeffsPerAxis) noexcept {
    auto ext = makeExtents(nCoeffsPerAxis);
    auto mapping = typename Mdspan::mapping_type{ext};
    return mapping.required_span_size();
}
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::buildCoeffs(int nCoeffsPerAxis, Func f,
                                                                    const DomainParams &dp) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto nodes = makeBuffer<Scalar, NCOEFFS>(nCoeffs);
    for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx)
        nodes[coeffIdx] =
            detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi / (2.0 * double(nCoeffsPerAxis)));

    std::array<int, DIM> extents{};
    extents.fill(nCoeffsPerAxis);

    // sample f on Chebyshev grid
    forEachIndex<DIM>(extents, [&](const std::array<int, DIM> &idx) {
        CanonicalInput domainPoint{};
        for (std::size_t d = 0; d < DIM; ++d) domainPoint[d] = nodes[static_cast<std::size_t>(idx[d])];
        const auto y = toCanonicalOutput(f(fromCanonicalInput(mapToDomain(domainPoint, dp))));
        for (std::size_t k = 0; k < OUT_DIM; ++k) coeff(idx, k) = y[k];
    });

    convertAxesToMonomialBasis(nCoeffsPerAxis, nodes);
    reverseCoefficientOrder(nCoeffsPerAxis);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::convertAxesToMonomialBasis(int nCoeffsPerAxis,
                                                                     const Buffer<Scalar, NCOEFFS> &nodes) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto rhs = makeBuffer<Scalar, NCOEFFS>(nCoeffs);
    auto alpha = makeBuffer<Scalar, NCOEFFS>(nCoeffs);
    auto mono = makeBuffer<Scalar, NCOEFFS>(nCoeffs);

    std::array<int, DIM> extents{};
    extents.fill(nCoeffsPerAxis);
    std::array<int, DIM> baseIndex{};
    for (std::size_t axis = 0; axis < DIM; ++axis) {
        auto innerExtents = extents;
        innerExtents[axis] = 1;
        forEachIndex<DIM>(innerExtents, [&](const std::array<int, DIM> &base) {
            for (std::size_t k = 0; k < OUT_DIM; ++k) {
                for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(coeffIdx);
                    rhs[coeffIdx] = coeff(baseIndex, k);
                }
                alpha = detail::bjorckPereyra<NCOEFFS, Scalar, Scalar>(nodes, rhs);
                mono = detail::newtonToMonomial<NCOEFFS, Scalar, Scalar>(alpha, nodes);
                for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(coeffIdx);
                    coeff(baseIndex, k) = mono[coeffIdx];
                }
            }
        });
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::reverseCoefficientOrder(int nCoeffsPerAxis) {
    std::array<int, DIM> extents{};
    extents.fill(nCoeffsPerAxis);
    std::array<int, DIM> baseIndex{};
    for (std::size_t axis = 0; axis < DIM; ++axis) {
        auto innerExtents = extents;
        innerExtents[axis] = 1;
        forEachIndex<DIM>(innerExtents, [&](const std::array<int, DIM> &base) {
            for (std::size_t k = 0; k < OUT_DIM; ++k) {
                int frontCoeff = 0;
                int backCoeff = nCoeffsPerAxis - 1;
                while (frontCoeff < backCoeff) {
                    baseIndex = base;
                    baseIndex[axis] = frontCoeff;
                    auto &a = coeff(baseIndex, k);
                    baseIndex[axis] = backCoeff;
                    auto &b = coeff(baseIndex, k);
                    std::swap(a, b);
                    ++frontCoeff;
                    --backCoeff;
                }
            }
        });
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<class Point>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::CanonicalInput FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::toCanonicalInput(
    const Point &x) noexcept {
    return detail::fixedContainerCast<CanonicalInput>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::InputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::fromCanonicalInput(const CanonicalInput &x) noexcept {
    return detail::fixedContainerCast<InputType>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::CanonicalOutput
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::toCanonicalOutput(const OutputType &x) noexcept {
    return detail::fixedContainerCast<CanonicalOutput>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::OutputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::fromCanonicalOutput(const CanonicalOutput &x) noexcept {
    return detail::fixedContainerCast<OutputType>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::CanonicalInput
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::mapToDomain(const CanonicalInput &x, const DomainParams &dp) noexcept {
    if (dp.identityDomain) return x;
    return polyfit::internal::helpers::mapToDomainArray<Scalar, DIM>(x, dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::CanonicalInput
FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::mapFromDomain(const CanonicalInput &x) const noexcept {
    if constexpr (FUSION_MODE == FusionMode::Always) {
        return x;
    } else {
        if (domain_.identityDomain) return x;
        return polyfit::internal::helpers::mapFromDomainArray<Scalar, DIM>(x, domain_.invSpan, domain_.sumEndpoints);
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::fuseNDDomain(DomainParams &dp, int nCoeffsPerAxis) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto fiber = makeBuffer<Scalar, NCOEFFS>(nCoeffs);
    std::array<int, DIM> extents{};
    extents.fill(nCoeffsPerAxis);
    std::array<int, DIM> baseIndex{};

    for (std::size_t axis = 0; axis < DIM; ++axis) {
        const auto alpha = Scalar(2) * static_cast<Scalar>(dp.invSpan[axis]);
        const auto beta = -static_cast<Scalar>(dp.sumEndpoints[axis]) * static_cast<Scalar>(dp.invSpan[axis]);
        // Axis already maps to [-1,1] — fuseLinearMap(alpha=1, beta=0) is
        // a numerical no-op; skip the ~n^2-per-fiber work.
        if (alpha == Scalar(1) && beta == Scalar(0)) {
            dp.invSpan[axis] = Scalar(0.5);
            dp.sumEndpoints[axis] = Scalar(0);
            continue;
        }

        auto innerExtents = extents;
        innerExtents[axis] = 1;
        forEachIndex<DIM>(innerExtents, [&](const std::array<int, DIM> &base) {
            for (std::size_t k = 0; k < OUT_DIM; ++k) {
                for (std::size_t i = 0; i < nCoeffs; ++i) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(i);
                    fiber[i] = coeff(baseIndex, k);
                }
                polyfit::internal::helpers::fuseLinearMap(fiber.data(), nCoeffs, alpha, beta);
                for (std::size_t i = 0; i < nCoeffs; ++i) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(i);
                    coeff(baseIndex, k) = fiber[i];
                }
            }
        });

        dp.invSpan[axis] = Scalar(0.5);
        dp.sumEndpoints[axis] = Scalar(0);
    }

    dp.identityDomain = polyfit::internal::helpers::isIdMap(dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::computeScaling(const InputType &a, const InputType &b,
                                                                      DomainParams &dp) const noexcept {
    polyfit::internal::helpers::computeScalingArray<Scalar, DIM>(toCanonicalInput(a), toCanonicalInput(b), dp.invSpan,
                                                                 dp.sumEndpoints);
    dp.identityDomain = polyfit::internal::helpers::isIdMap(dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE, ScalarKernel SK, std::size_t HK>
template<std::size_t Rank, class F>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE, SK, HK>::forEachIndex(const std::array<int, Rank> &ext, F &&body) {
    std::array<int, Rank> idx{};
    while (true) {
        body(idx);
        for (std::size_t d = 0; d < Rank; ++d) {
            if (++idx[d] < ext[d]) break;
            if (d == Rank - 1) return;
            idx[d] = 0;
        }
    }
}
} // namespace poly_eval

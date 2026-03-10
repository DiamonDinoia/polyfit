#pragma once

#include <cassert>
#include <stdexcept>
#include <xsimd/xsimd.hpp>

#include "helpers.h"
#include "macros.h"
#include "poly_eval.h"
#include "utils.h"

#include <poet/poet.hpp>

namespace poly_eval {
namespace detail {

constexpr int validate_positive_nCoeffs(const int n) {
    if (n <= 0) {
        throw std::invalid_argument("Runtime coefficient count must be positive (nCoeffs > 0)");
    }
    return n;
}

template<typename T> constexpr void validate_domain(const T &a, const T &b) {
    if constexpr (has_tuple_size_v<T>) {
        for (std::size_t i = 0; i < std::tuple_size_v<T>; ++i) {
            if (a[i] == b[i]) {
                throw std::invalid_argument("Domain endpoints must differ in every dimension");
            }
        }
    } else {
        if (a == b) {
            throw std::invalid_argument("Domain endpoints must differ");
        }
    }
}

} // namespace detail

// -----------------------------------------------------------------------------
// FuncEval Implementation (Runtime)
// -----------------------------------------------------------------------------

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<std::size_t CurrentNCoeffs, typename>
PF_C20CONSTEXPR FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::FuncEval(Func F, const int nCoeffs,
                                                                                      const InputType a,
                                                                                      const InputType b,
                                                                                      const InputType *pts)
    : invSpan(InputType(1) / (b - a)), sumEndpoints(b + a) {
    detail::validate_domain(a, b);
    const auto validatedNCoeffs = detail::validate_positive_nCoeffs(nCoeffs);
    coeffsBuf.resize(static_cast<std::size_t>(validatedNCoeffs));
    initializeCoeffs(F, pts);
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<std::size_t CurrentNCoeffs, typename>
PF_C20CONSTEXPR FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::FuncEval(Func F, const InputType a,
                                                                                      const InputType b,
                                                                                      const InputType *pts)
    : invSpan(InputType(1) / (b - a)), sumEndpoints(b + a) {
    static_assert(CurrentNCoeffs > 0, "Compile-time coefficient count must be positive");
    detail::validate_domain(a, b);
    initializeCoeffs(F, pts);
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<bool>
constexpr typename FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::OutputType PF_ALWAYS_INLINE
FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::operator()(const InputType pt) const noexcept {
    const auto xi = mapFromDomain(pt);
    return horner<NCoeffsCt>(xi, coeffsBuf.data(), coeffsBuf.size()); // Pass data pointer and size
}

// Batch evaluation implementation using SIMD and unrolling
PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<int OuterUnrollFactor, bool ptsAligned, bool outAligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::evalBatch(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
    return horner<NCoeffsCt, ptsAligned, outAligned, OuterUnrollFactor>(
        pts, out, numPoints, coeffsBuf.data(), coeffsBuf.size(),
        [this](const auto v) { return this->mapFromDomain(v); });
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<bool ptsAligned, bool outAligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::operator()(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
#ifdef PF_OUTER_UNROLL
    PF_C23STATIC constexpr auto unrollFactor = PF_OUTER_UNROLL;
#else
    PF_C23STATIC constexpr auto unrollFactor = 0;
#endif

    if constexpr (ptsAligned) {
        if constexpr (outAligned) {
            return evalBatch<unrollFactor, true, true>(pts, out, numPoints);
        } else {
            return evalBatch<unrollFactor, true, false>(pts, out, numPoints);
        }
    }

    return evalBatch<unrollFactor, false, false>(pts, out, numPoints);
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
PF_C20CONSTEXPR const Buffer<typename FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::OutputType,
                             NCoeffsCt> &
FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::coeffs() const noexcept {
    return coeffsBuf;
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
constexpr std::size_t FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::nCoeffs() const noexcept {
    return coeffsBuf.size();
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<class T>
[[nodiscard]] PF_ALWAYS_INLINE constexpr T
FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::mapToDomain(const T value) const noexcept {
    return polyfit::internal::helpers::map_to_domain_scalar(value, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
template<class T>
[[nodiscard]] PF_ALWAYS_INLINE constexpr T
FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::mapFromDomain(const T value) const noexcept {
    if constexpr (Fusion == FusionMode::always)
        return value;
    else if constexpr (Fusion == FusionMode::never)
        return polyfit::internal::helpers::map_from_domain_scalar(value, invSpan, sumEndpoints);
    else {
        if (domainFused) return value;
        return polyfit::internal::helpers::map_from_domain_scalar(value, invSpan, sumEndpoints);
    }
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
PF_C20CONSTEXPR void FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::initializeCoeffs(
    Func F, const InputType *pts) {
    // 1) allocate
    auto grid = make_buffer<InputType, NCoeffsCt>(coeffsBuf.size());
    auto samples = make_buffer<OutputType, NCoeffsCt>(coeffsBuf.size());

    // 2) fill
    const auto nCoeffs = coeffsBuf.size();
    for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
        grid[coeffIdx] = pts ? pts[coeffIdx]
                             : InputType(detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi /
                                                     (2.0 * double(nCoeffs))));
    }
    for (std::size_t sampleIdx = 0; sampleIdx < nCoeffs; ++sampleIdx) {
        samples[sampleIdx] = F(mapToDomain(grid[sampleIdx]));
    }

    // 3) compute Newton → monomial
    auto newton = detail::bjorck_pereyra<NCoeffsCt, InputType, OutputType>(grid, samples);
    auto temp_monomial = detail::newton_to_monomial<NCoeffsCt, InputType, OutputType>(newton, grid);
    assert(temp_monomial.size() == coeffsBuf.size() && "size mismatch!");

    std::copy(temp_monomial.begin(), temp_monomial.end(), coeffsBuf.begin());

    // 4) optional refine
    refine(grid, samples);

    // 5) Domain fusion: fuse the linear domain mapping into polynomial coefficients.
    //    p(t) → q(x) = p(alpha*x + beta) where alpha = 2*invSpan, beta = -sumEndpoints*invSpan.
    //    After this, evaluation skips per-point mapping.
    if constexpr (Fusion != FusionMode::never) {
#if __cplusplus >= 202602L
        {
#else
        PF_IF_NOT_CONSTEVAL {
#endif
            using Scalar = typename value_type_or_identity<InputType>::type;
            const auto alpha = Scalar(2) * static_cast<Scalar>(invSpan);
            const auto beta = -static_cast<Scalar>(sumEndpoints) * static_cast<Scalar>(invSpan);
            const auto coeffCountForFusion = static_cast<int>(coeffsBuf.size());

            bool should_fuse;
            if constexpr (Fusion == FusionMode::always) {
                should_fuse = true;
            } else { // FusionMode::auto_
                const auto cond_base = detail::math::abs(alpha) + detail::math::abs(beta) + Scalar(1);
                constexpr auto max_log = Scalar(std::numeric_limits<Scalar>::digits10 - 3);
                should_fuse = (coeffCountForFusion > 1 &&
                               Scalar(coeffCountForFusion - 1) * detail::math::log10(cond_base) < max_log);
            }

            if (should_fuse) {
                polyfit::internal::helpers::fuse_linear_map(coeffsBuf.data(), coeffsBuf.size(), alpha, beta);
                invSpan = InputType(0.5);
                sumEndpoints = InputType(0);
                if constexpr (Fusion == FusionMode::auto_) domainFused = true;
            }
        }
    }
}

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
PF_C20CONSTEXPR void
FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::refine(const Buffer<InputType, NCoeffsCt> &chebNodes,
                                                                    const Buffer<OutputType, NCoeffsCt> &samples) {
    const auto nCoeffs = coeffsBuf.size();
    std::reverse(coeffsBuf.begin(), coeffsBuf.end()); // to Horner order once

    // Compensated Horner gives more accurate residuals, preventing refinement
    // divergence at high coefficient counts where standard Horner rounding errors exceed
    // the correction magnitude from bjorck_pereyra + newton_to_monomial.
    // Extra passes for n > 32: the compensated residuals converge where standard ones diverge.
    const std::size_t total_iters = RefineIters + (nCoeffs > 32 ? 2 : 0);

    for (std::size_t pass = 0; pass < total_iters; ++pass) {
        auto residuals = make_buffer<OutputType, NCoeffsCt>(nCoeffs);
        for (std::size_t sampleIdx = 0; sampleIdx < nCoeffs; ++sampleIdx) {
            auto p_val = poly_eval::compensated_horner<NCoeffsCt>(chebNodes[sampleIdx], coeffsBuf.data(), nCoeffs);
            residuals[sampleIdx] = samples[sampleIdx] - p_val;
        }
        auto newton_r = detail::bjorck_pereyra<NCoeffsCt, InputType, OutputType>(chebNodes, residuals);
        auto mono_r = detail::newton_to_monomial<NCoeffsCt, InputType, OutputType>(newton_r, chebNodes);

        // mono_r is low-order-first; coeffsBuf is highest-order-first — add reversed
        for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
            coeffsBuf[nCoeffs - 1 - coeffIdx] += mono_r[coeffIdx];
        }
    }
    // coeffsBuf are now in Horner order (highest-order-first) — done
}

// ------------------------------ truncate -------------------------------

template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
PF_C20CONSTEXPR void FuncEval<Func, NCoeffsCt, RefineIters, Fusion>::truncate(
    typename value_type_or_identity<OutputType>::type eps) noexcept {
    if constexpr (NCoeffsCt == 0) {
        // coeffsBuf are in Horner order: [0]=highest-order term, [size-1]=constant
        std::size_t skip = 0;
        while (skip + 1 < coeffsBuf.size() && std::abs(coeffsBuf[skip]) < eps) ++skip;
        if (skip > 0) coeffsBuf.erase(coeffsBuf.begin(), coeffsBuf.begin() + static_cast<std::ptrdiff_t>(skip));
    }
}

// ------------------------------ ctor -----------------------------------

template<typename... EvalTypes> PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const EvalTypes &...evals) {
    /* Copy per‑poly scaling */
    invSpan = {evals.invSpan...};
    sumEndpoints = {evals.sumEndpoints...};

    /* Coefficient-count-dependent storage */
    if constexpr (maxNCoeffsCt == 0) {
        const std::size_t maxNCoeffs = std::max({std::size_t(evals.coeffsBuf.size())...});
        coeffStore.assign(kFPad * maxNCoeffs, OutputType{});
        coeffs = decltype(coeffs){coeffStore.data(), maxNCoeffs, kFPad};
    } else {
        coeffs = decltype(coeffs){coeffStore.data(), maxNCoeffsCt, kFPad};
    }

    /* Copy coefficients and pad */
    copyCoeffs<0>(evals...);
    zeroPadCoeffs();
}

template<typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const FuncEvalMany &other)
    : coeffStore(other.coeffStore), coeffs{coeffStore.data(), other.coeffs.extents()}, invSpan(other.invSpan),
      sumEndpoints(other.sumEndpoints) {}

template<typename... EvalTypes>
PF_C20CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(const FuncEvalMany &other) -> FuncEvalMany & {
    if (this != &other) {
        coeffStore = other.coeffStore;
        invSpan = other.invSpan;
        sumEndpoints = other.sumEndpoints;
        coeffs = decltype(coeffs){coeffStore.data(), other.coeffs.extents()};
    }
    return *this;
}

template<typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(FuncEvalMany &&other) noexcept
    : coeffStore(std::move(other.coeffStore)), coeffs{coeffStore.data(), other.coeffs.extents()},
      invSpan(std::move(other.invSpan)), sumEndpoints(std::move(other.sumEndpoints)) {}

template<typename... EvalTypes>
PF_C20CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(FuncEvalMany &&other) noexcept -> FuncEvalMany & {
    if (this != &other) {
        coeffStore = std::move(other.coeffStore);
        invSpan = std::move(other.invSpan);
        sumEndpoints = std::move(other.sumEndpoints);
        coeffs = decltype(coeffs){coeffStore.data(), other.coeffs.extents()};
    }
    return *this;
}

// ------------------------------ size / coefficient count --------------------------

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::size() const noexcept { return kF; }

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::nCoeffs() const noexcept {
    return coeffs.extent(0);
}

// ------------------------------ scalar eval ----------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(InputType x) const noexcept -> std::array<OutputType, kF> {
    alignas(kAlignment) std::array<InputType, kFPad> xu{};
    poet::static_for<kF>([&](auto i) { xu[i] = xsimd::fms(InputType(2.0), x, sumEndpoints[i]) * invSpan[i]; });

    alignas(kAlignment) std::array<OutputType, kFPad> res{};
    horner_transposed<kFPad, maxNCoeffsCt, vectorWidth, true>(
        xu.data(), coeffs.data_handle(), res.data(), kFPad, static_cast<std::size_t>(coeffs.extent(0)));
    return extractReal(res);
}
PF_FAST_EVAL_END

// ------------------------------ array eval -----------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(const std::array<InputType, kF> &xs) const noexcept
    -> std::array<OutputType, kF> {
    alignas(kAlignment) std::array<InputType, kFPad> xu{};
    poet::static_for<kF>([&](auto i) { xu[i] = xsimd::fms(InputType(2.0), xs[i], sumEndpoints[i]) * invSpan[i]; });

    alignas(kAlignment) std::array<OutputType, kFPad> res{};
    horner_transposed<kFPad, maxNCoeffsCt, vectorWidth, true>(
        xu.data(), coeffs.data_handle(), res.data(), kFPad, static_cast<std::size_t>(coeffs.extent(0)));
    return extractReal(res);
}
PF_FAST_EVAL_END

// ------------------------------ bulk eval ------------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::operator()(const InputType *PF_RESTRICT x, OutputType *PF_RESTRICT out,
                                            std::size_t numPoints) const noexcept {
    PF_C23STATIC constexpr std::size_t M = kF;
    PF_C23STATIC constexpr std::size_t simdSize = xsimd::batch<InputType>::size;
    PF_C23STATIC constexpr std::size_t stride = kFPad;

    // UF independent Horner chains for ILP, tuned to vector register pressure.
    // Explicit multi-accumulator: shared coeff broadcast across UF lanes per k-step.
    PF_C23STATIC constexpr std::size_t UF = detail::optimal_many_eval_uf<OutputType>();
    PF_C23STATIC constexpr std::size_t block = simdSize * UF;

    const auto nCoeffs = static_cast<std::size_t>(coeffs.extent(0));

    // Horner coefficient loop: CT → static_for (full unrolling), RT → plain for.
    auto forEachCoeff = [&](auto step) {
        if constexpr (maxNCoeffsCt != 0)
            poet::static_for<1, maxNCoeffsCt>(step);
        else
            for (std::size_t k = 1; k < nCoeffs; ++k) step(k);
    };

    poet::static_for<M>([&](auto m) {
        const auto invSpanValue = invSpan[m];
        const auto sumEndpointsValue = sumEndpoints[m];
        const OutputType *col_m = coeffs.data_handle() + std::size_t(m);

        auto scatter = [&](auto acc, std::size_t base) {
            alignas(xsimd::batch<OutputType>::arch_type::alignment()) OutputType tmp[simdSize];
            acc.store_aligned(tmp);
            poet::static_for<simdSize>([&](auto s) {
                out[(base + std::size_t(s)) * M + m] = tmp[s];
            });
        };

        const auto twoVec = xsimd::batch<InputType>(InputType(2.0));
        const auto sumEndpointsVec = xsimd::batch<InputType>(sumEndpointsValue);
        const auto invSpanVec = xsimd::batch<InputType>(invSpanValue);
        auto mapSimd = [&](auto xv) { return xsimd::fms(twoVec, xv, sumEndpointsVec) * invSpanVec; };

        const auto tileEnd = detail::round_down<block>(numPoints);
        const auto simdEnd = detail::round_down<simdSize>(numPoints);

        for (std::size_t i = 0; i < tileEnd; i += block) {
            xsimd::batch<InputType> pt[UF];
            xsimd::batch<OutputType> acc[UF];

            poet::static_for<UF>([&](auto j) {
                pt[j] = mapSimd(xsimd::load_unaligned(x + i + j * simdSize));
                acc[j] = xsimd::batch<OutputType>(col_m[0]);
            });

            forEachCoeff([&](auto k) {
                const auto ck = xsimd::batch<OutputType>(col_m[std::size_t(k) * stride]);
                poet::static_for<UF>([&](auto j) { acc[j] = detail::fma(acc[j], pt[j], ck); });
            });

            poet::static_for<UF>([&](auto j) { scatter(acc[j], i + j * simdSize); });
        }

        for (std::size_t i = tileEnd; i < simdEnd; i += simdSize) {
            auto xu = mapSimd(xsimd::load_unaligned(x + i));
            auto acc = xsimd::batch<OutputType>(col_m[0]);

            forEachCoeff([&](auto k) {
                acc = detail::fma(acc, xu, xsimd::batch<OutputType>(col_m[std::size_t(k) * stride]));
            });

            scatter(acc, i);
        }

        for (std::size_t i = simdEnd; i < numPoints; ++i) {
            auto xu = (InputType(2.0) * x[i] - sumEndpointsValue) * invSpanValue;
            OutputType acc = col_m[0];

            forEachCoeff([&](auto k) {
                acc = detail::fma(acc, xu, col_m[std::size_t(k) * stride]);
            });

            out[i * M + m] = acc;
        }
    });
}
PF_FAST_EVAL_END

// ------------------------------ variadic convenience -------------------

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(InputType first, Ts... rest) const noexcept -> std::array<OutputType, kF> {
    static_assert(sizeof...(Ts) + 1 == kF, "Incorrect number of arguments");
    return operator()(std::array<InputType, kF>{first, static_cast<InputType>(rest)...});
}

// ------------------------------ tuple convenience ----------------------

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(const std::tuple<Ts...> &tup) const noexcept -> std::array<OutputType, kF> {
    static_assert(sizeof...(Ts) == kF, "Tuple size must equal number of polynomials");
    std::array<InputType, kF> xs{};
    std::apply([&](auto &&...e) { xs = {static_cast<InputType>(e)...}; }, tup);
    return operator()(xs);
}

// ------------------------------ copyCoeffs ----------------------------

template<typename... EvalTypes>
template<std::size_t I, typename FE, typename... Rest>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::copyCoeffs(const FE &eval, const Rest &...rest) {
    for (std::size_t k = 0; k < eval.coeffsBuf.size(); ++k) coeffs(k, I) = eval.coeffsBuf[k];
    for (std::size_t k = eval.coeffsBuf.size(); k < coeffs.extent(0); ++k) coeffs(k, I) = OutputType{};
    if constexpr (I + 1 < kF) copyCoeffs<I + 1>(rest...);
}

// ------------------------------ zeroPadCoeffs ------------------------

template<typename... EvalTypes> PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::zeroPadCoeffs() {
    for (std::size_t j = kF; j < kFPad; ++j)
        for (std::size_t k = 0; k < coeffs.extent(0); ++k) coeffs(k, j) = OutputType{};
}

// ------------------------------ truncate -------------------------------

template<typename... EvalTypes>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::truncate(typename value_type_or_identity<OutputType>::type eps) {
    // Scan from the leading coefficient row downward.
    std::size_t newNCoeffs = coeffs.extent(0);
    while (newNCoeffs > 1) {
        OutputType rowMax{};
        for (std::size_t j = 0; j < kF; ++j) rowMax = std::max(rowMax, std::abs(coeffs(newNCoeffs - 1, j)));
        if (rowMax >= eps) break;
        --newNCoeffs;
    }
    if (newNCoeffs < coeffs.extent(0)) {
        coeffs = decltype(coeffs){coeffStore.data(), newNCoeffs, kFPad};
    }
}

// ------------------------------ extractReal ---------------------------

template<typename... EvalTypes>
constexpr auto FuncEvalMany<EvalTypes...>::extractReal(const std::array<OutputType, kFPad> &full) const noexcept
    -> std::array<OutputType, kF> {
    if constexpr (kF == kFPad) {
        return full; // no padding
    }
    std::array<OutputType, kF> out{};
    poet::static_for<kF>([&](auto i) { out[i] = full[i]; });
    return out;
}

// Constructor for static coefficient count
template<class Func, std::size_t NCoeffsCt>
template<std::size_t C, typename>
constexpr FuncEvalND<Func, NCoeffsCt>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : coeffsFlat(), coeffsMd{coeffsFlat.data(), Extents{}} {
    detail::validate_domain(a, b);
    computeScaling(a, b);
    initialize(static_cast<int>(NCoeffsCt), f);
}

// Constructor for dynamic coefficient count
template<class Func, std::size_t NCoeffsCt>
template<std::size_t C, typename>
constexpr FuncEvalND<Func, NCoeffsCt>::FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b)
    : coeffsFlat(storage_required(detail::validate_positive_nCoeffs(nCoeffsPerAxis))),
      coeffsMd{coeffsFlat.data(), makeExtents(nCoeffsPerAxis)} {
    detail::validate_domain(a, b);
    computeScaling(a, b);
    initialize(nCoeffsPerAxis, f);
}

template<class Func, std::size_t NCoeffsCt>
constexpr std::size_t FuncEvalND<Func, NCoeffsCt>::nCoeffsPerAxis() const noexcept {
    return static_cast<std::size_t>(coeffsMd.extent(0));
}

template<class Func, std::size_t NCoeffsCt>
FuncEvalND<Func, NCoeffsCt>::FuncEvalND(const FuncEvalND &other)
    : invSpan(other.invSpan), sumEndpoints(other.sumEndpoints), coeffsFlat(other.coeffsFlat),
      coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCoeffsCt>
auto FuncEvalND<Func, NCoeffsCt>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        invSpan = other.invSpan;
        sumEndpoints = other.sumEndpoints;
        coeffsFlat = other.coeffsFlat;
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

template<class Func, std::size_t NCoeffsCt>
FuncEvalND<Func, NCoeffsCt>::FuncEvalND(FuncEvalND &&other) noexcept
    : invSpan(std::move(other.invSpan)), sumEndpoints(std::move(other.sumEndpoints)),
      coeffsFlat(std::move(other.coeffsFlat)), coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCoeffsCt>
auto FuncEvalND<Func, NCoeffsCt>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        invSpan = std::move(other.invSpan);
        sumEndpoints = std::move(other.sumEndpoints);
        coeffsFlat = std::move(other.coeffsFlat);
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

// Evaluate via Horner's method
PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCoeffsCt>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCoeffsCt>::OutputType
FuncEvalND<Func, NCoeffsCt>::operator()(const InputType &x) const {
    const int nCoeffsRt = (NCoeffsCt ? static_cast<int>(NCoeffsCt) : static_cast<int>(coeffsMd.extent(0)));
    return poly_eval::horner<NCoeffsCt, SIMD, OutputType>(mapFromDomain(x), coeffsMd, nCoeffsRt);
}
PF_FAST_EVAL_END

// coeff_impl
template<class Func, std::size_t NCoeffsCt>
template<typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, NCoeffsCt>::Scalar &
FuncEvalND<Func, NCoeffsCt>::coeffImpl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffsMd(static_cast<std::size_t>(idx[I])..., k);
}

template<class Func, std::size_t NCoeffsCt>
template<class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCoeffsCt>::Scalar &
FuncEvalND<Func, NCoeffsCt>::coeff(const IdxArray &idx, std::size_t k) noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<dim>{});
}

template<class Func, std::size_t NCoeffsCt>
auto FuncEvalND<Func, NCoeffsCt>::makeExtents(int nCoeffsPerAxis) noexcept -> Extents {
    if constexpr (isStatic) {
        return detail::make_static_extents<NCoeffsCt, dim, outDim>(std::make_index_sequence<dim>{});
    } else {
        return makeExtents(nCoeffsPerAxis, std::make_index_sequence<dim + 1>{});
    }
}

template<class Func, std::size_t NCoeffsCt>
template<std::size_t... Is>
auto FuncEvalND<Func, NCoeffsCt>::makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept -> Extents {
    return Extents{(Is < dim ? static_cast<std::size_t>(nCoeffsPerAxis) : static_cast<std::size_t>(outDim))...};
}

template<class Func, std::size_t NCoeffsCt>
constexpr std::size_t FuncEvalND<Func, NCoeffsCt>::storage_required(const int nCoeffsPerAxis) noexcept {
    auto ext = makeExtents(nCoeffsPerAxis);
    auto mapping = typename Mdspan::mapping_type{ext};
    return mapping.required_span_size();
}
template<class Func, std::size_t NCoeffsCt>
constexpr void FuncEvalND<Func, NCoeffsCt>::initialize(int nCoeffsPerAxis, Func f) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto nodes = make_buffer<Scalar, NCoeffsCt>(nCoeffs);
    for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx)
        nodes[coeffIdx] =
            detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi / (2.0 * double(nCoeffsPerAxis)));

    std::array<int, dim> extents{};
    extents.fill(nCoeffsPerAxis);

    // sample f on Chebyshev grid
    forEachIndex<dim>(extents, [&](const std::array<int, dim> &idx) {
        InputType domainPoint{};
        for (std::size_t d = 0; d < dim; ++d) domainPoint[d] = nodes[static_cast<std::size_t>(idx[d])];
        OutputType y = f(mapToDomain(domainPoint));
        for (std::size_t k = 0; k < outDim; ++k) coeff(idx, k) = y[k];
    });

    convertAxesToMonomial(nCoeffsPerAxis, nodes);
    reverseAxes(nCoeffsPerAxis);
}

template<class Func, std::size_t NCoeffsCt>
constexpr void FuncEvalND<Func, NCoeffsCt>::convertAxesToMonomial(int nCoeffsPerAxis,
                                                                            const Buffer<Scalar, NCoeffsCt> &nodes) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto rhs = make_buffer<Scalar, NCoeffsCt>(nCoeffs);
    auto alpha = make_buffer<Scalar, NCoeffsCt>(nCoeffs);
    auto mono = make_buffer<Scalar, NCoeffsCt>(nCoeffs);

    std::array<int, dim> extents{};
    extents.fill(nCoeffsPerAxis);
    std::array<int, dim> baseIndex{};
    for (std::size_t axis = 0; axis < dim; ++axis) {
        auto innerExtents = extents;
        innerExtents[axis] = 1;
        forEachIndex<dim>(innerExtents, [&](const std::array<int, dim> &base) {
            for (std::size_t k = 0; k < outDim; ++k) {
                for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(coeffIdx);
                    rhs[coeffIdx] = coeff(baseIndex, k);
                }
                alpha = detail::bjorck_pereyra<NCoeffsCt, Scalar, Scalar>(nodes, rhs);
                mono = detail::newton_to_monomial<NCoeffsCt, Scalar, Scalar>(alpha, nodes);
                for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
                    baseIndex = base;
                    baseIndex[axis] = static_cast<int>(coeffIdx);
                    coeff(baseIndex, k) = mono[coeffIdx];
                }
            }
        });
    }
}

template<class Func, std::size_t NCoeffsCt>
constexpr void FuncEvalND<Func, NCoeffsCt>::reverseAxes(int nCoeffsPerAxis) {
    std::array<int, dim> extents{};
    extents.fill(nCoeffsPerAxis);
    std::array<int, dim> baseIndex{};
    for (std::size_t axis = 0; axis < dim; ++axis) {
        auto innerExtents = extents;
        innerExtents[axis] = 1;
        forEachIndex<dim>(innerExtents, [&](const std::array<int, dim> &base) {
            for (std::size_t k = 0; k < outDim; ++k) {
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

template<class Func, std::size_t NCoeffsCt>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCoeffsCt>::InputType
FuncEvalND<Func, NCoeffsCt>::mapToDomain(const InputType &x) const noexcept {
    return polyfit::internal::helpers::map_to_domain_array<Scalar, dim>(x, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCoeffsCt>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCoeffsCt>::InputType
FuncEvalND<Func, NCoeffsCt>::mapFromDomain(const InputType &x) const noexcept {
    return polyfit::internal::helpers::map_from_domain_array<Scalar, dim>(x, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCoeffsCt>
constexpr void FuncEvalND<Func, NCoeffsCt>::computeScaling(const InputType &a, const InputType &b) noexcept {
    polyfit::internal::helpers::compute_scaling_array<Scalar, dim>(a, b, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCoeffsCt>
// Odometer-style multi-index iteration: visits every index tuple in the
// Cartesian product [0, ext[0]) × [0, ext[1]) × … × [0, ext[Rank-1]).
// The innermost (d=0) dimension increments fastest, like a little-endian
// counter. Each visited index tuple is passed to `body`.
template<std::size_t Rank, class F>
constexpr void FuncEvalND<Func, NCoeffsCt>::forEachIndex(const std::array<int, Rank> &ext, F &&body) {
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

// -----------------------------------------------------------------------------
// make_func_eval API implementations (tag-based, C++17+)
// -----------------------------------------------------------------------------

// Compile-time coefficient count (1D or ND)
template<std::size_t NCoeffsCt, class Func, class... Tags,
         std::enable_if_t<(NCoeffsCt > 0) && detail::all_tags_v<Tags...>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...) {
    using InputType = typename function_traits<Func>::arg0_type;
    constexpr std::size_t refineIters = detail::extract_tag<iters, 1, Tags...>();
    constexpr FusionMode fusionMode = detail::extract_fusion_v<Tags...>;
    if constexpr (has_tuple_size_v<std::remove_cvref_t<InputType>>) {
        return FuncEvalND<Func, NCoeffsCt>(F, a, b);
    } else {
        return FuncEval<Func, NCoeffsCt, refineIters, fusionMode>(F, a, b);
    }
}

// Runtime coefficient count (1D or ND)
template<class Func, typename IntType, class... Tags,
         std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>> && detail::all_tags_v<Tags...>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto
make_func_eval(Func F, IntType nCoeffs, typename function_traits<Func>::arg0_type a,
               typename function_traits<Func>::arg0_type b, Tags...) {
    using InputType = typename function_traits<Func>::arg0_type;
    using RawInputType = std::remove_cvref_t<InputType>;
    constexpr std::size_t refineIters = detail::extract_tag<iters, 1, Tags...>();
    constexpr FusionMode fusionMode = detail::extract_fusion_v<Tags...>;

    if constexpr (has_tuple_size_v<RawInputType>) {
        return FuncEvalND<Func, 0>(F, detail::validate_positive_nCoeffs(static_cast<int>(nCoeffs)), a, b);
    } else {
        return FuncEval<Func, 0, refineIters, fusionMode>(F, detail::validate_positive_nCoeffs(static_cast<int>(nCoeffs)),
                                                          a, b);
    }
}

// Runtime error tolerance (1D or ND) — search for the first coefficient count that meets eps
template<class Func, typename FloatType, class... Tags,
         std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>> && detail::all_tags_v<Tags...>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...) {
    using RawInputType = std::remove_cvref_t<typename function_traits<Func>::arg0_type>;
    constexpr std::size_t maxNCoeffsVal = detail::extract_tag<maxNCoeffs, 32, Tags...>();
    constexpr std::size_t numEvalPointsVal = detail::extract_tag<evalPts, 100, Tags...>();
    constexpr std::size_t refineIters = detail::extract_tag<iters, 1, Tags...>();
    constexpr FusionMode fusionMode = detail::extract_fusion_v<Tags...>;

    using Evaluator =
        std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, 0>, FuncEval<Func, 0, refineIters, fusionMode>>;
    const auto evalPoints = detail::linspace(a, b, int(numEvalPointsVal));

    auto computeMaxErr = [&](const Evaluator &evaluator) {
        double maxErr = 0.0;
        for (const auto &pt : evalPoints) {
            maxErr = std::max(detail::relative_l2_norm(F(pt), evaluator(pt)), maxErr);
        }
        return maxErr;
    };

    if (eps <= FloatType(0)) {
        throw std::invalid_argument("Requested error tolerance must be positive");
    }

    detail::validate_domain(a, b);

    for (std::size_t candidateNCoeffs = 1; candidateNCoeffs <= maxNCoeffsVal; ++candidateNCoeffs) {
        Evaluator evaluator(F, int(candidateNCoeffs), a, b);
        if (computeMaxErr(evaluator) <= eps) {
            return evaluator;
        }
    }

    const Evaluator maxNCoeffsEvaluator(F, int(maxNCoeffsVal), a, b);
    throw std::runtime_error("No coefficient count found for requested error tolerance. eps=" + std::to_string(eps) +
                             ", maxNCoeffs=" + std::to_string(maxNCoeffsVal) +
                             ", maxErr=" + std::to_string(computeMaxErr(maxNCoeffsEvaluator)));
}

// Function pointer overload
template<std::size_t NCoeffsCt, class... Tags, typename Func,
         std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>> && detail::all_tags_v<Tags...>,
                          int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                                  typename function_traits<Func *>::arg0_type b, Tags...) {
    using InputType = typename function_traits<Func *>::arg0_type;
    constexpr std::size_t refineIters = detail::extract_tag<iters, 1, Tags...>();
    constexpr FusionMode fusionMode = detail::extract_fusion_v<Tags...>;
    if constexpr (has_tuple_size_v<std::remove_cvref_t<InputType>>) {
        auto funcWrapper = [f](const InputType &in) { return f(in); };
        return FuncEvalND<decltype(funcWrapper), NCoeffsCt>(funcWrapper, a, b);
    } else {
        return FuncEval<Func *, NCoeffsCt, refineIters, fusionMode>(f, a, b, nullptr);
    }
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double epsVal, auto a, auto b, std::size_t maxNCoeffsVal, std::size_t numEvalPointsVal,
         std::size_t refineIters, class Func>
[[nodiscard]] constexpr auto make_func_eval(Func F) {
    using RawInputType = std::remove_cvref_t<typename function_traits<Func>::arg0_type>;
    static_assert(maxNCoeffsVal > 0, "Max coefficient count must be positive.");
    static_assert(numEvalPointsVal > 1, "Number of evaluation points must be greater than 1.");

    // Recursive constexpr search for the minimal coefficient count.
    constexpr auto nCoeffs = [F] {
        constexpr auto computeError = [F](const auto &evaluator) {
            constexpr auto ep = detail::linspace<static_cast<int>(numEvalPointsVal)>(a, b);
            double maxErr = 0.0;
            for (const auto &pt : ep) {
                const auto actual = F(pt);
                const auto approx = evaluator.template operator()<false>(pt);
                maxErr = std::max(detail::relative_l2_norm(actual, approx), maxErr);
            }
            return maxErr;
        };
        int result = 0;
        poet::static_for<1, maxNCoeffsVal + 1>([&](auto i) {
            if (result != 0) return;
            using Evaluator = std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, i>,
                                                 FuncEval<Func, i, refineIters>>;
            if constexpr (computeError(Evaluator(F, a, b)) <= epsVal) {
                result = i;
            }
        });
        return result;
    }();
    static_assert(nCoeffs != 0, "No coefficient count found for requested error tolerance.");
    using Evaluator = std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, nCoeffs>,
                                         FuncEval<Func, nCoeffs, refineIters>>;
    return Evaluator(F, a, b);
}
#endif // PF_HAS_CONSTEXPR_EPS_OVERLOAD

template<typename... EvalTypes>
[[nodiscard]] PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval

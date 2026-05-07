#pragma once

#include "impl_common.hpp"

// Out-of-class definitions for poly_eval::FuncEval.

namespace poly_eval {

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR FuncEval<Func, NCOEFFS, ITERS, FUSION>::FuncEval(Func F, const int nCoeffs, const InputType a,
                                                                 const InputType b, const InputType *pts) {
    DomainParams dp;
    dp.invSpan = InputType(1) / (b - a);
    dp.sumEndpoints = b + a;
    dp.identityMap = polyfit::internal::helpers::isIdMap(dp.invSpan, dp.sumEndpoints);
    initialize(detail::RuntimeCountTag{}, F, nCoeffs, a, b, pts, dp);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR FuncEval<Func, NCOEFFS, ITERS, FUSION>::FuncEval(Func F, const InputType a, const InputType b,
                                                                 const InputType *pts) {
    DomainParams dp;
    dp.invSpan = InputType(1) / (b - a);
    dp.sumEndpoints = b + a;
    dp.identityMap = polyfit::internal::helpers::isIdMap(dp.invSpan, dp.sumEndpoints);
    initialize(detail::CompileTimeCountTag{}, F, a, b, pts, dp);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void
FuncEval<Func, NCOEFFS, ITERS, FUSION>::initialize(detail::CompileTimeCountTag, Func F, const InputType a,
                                                   const InputType b, const InputType *pts, DomainParams &dp) {
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    detail::validateDomain(a, b);
    initializeCoeffs(F, pts, dp);
    if constexpr (kStoresDomain) domain_ = dp;
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::initialize(detail::RuntimeCountTag, Func F,
                                                                         const int nCoeffs, const InputType a,
                                                                         const InputType b, const InputType *pts,
                                                                         DomainParams &dp) {
    static_assert(NCOEFFS == 0, "Runtime coefficient count is only valid for runtime-sized evaluators");
    detail::validateDomain(a, b);
    const auto validatedNCoeffs = detail::validatePositiveCoeffCount(nCoeffs);
    coeffsBuf.resize(static_cast<std::size_t>(validatedNCoeffs));
    initializeCoeffs(F, pts, dp);
    if constexpr (kStoresDomain) domain_ = dp;
}

// ---------- ND runtime init: same one-shot fusion gating as compile-time ----

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool>
constexpr typename FuncEval<Func, NCOEFFS, ITERS, FUSION>::OutputType PF_ALWAYS_INLINE
FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(const InputType pt) const noexcept {
    const auto xi = mapFromDomain(pt);
    return hybrid<NCOEFFS>(xi, coeffsBuf.data(), coeffsBuf.size());
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool Dummy, class V>
constexpr auto PF_ALWAYS_INLINE FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(V pt) const noexcept
    -> enable_if_t<!std::is_same_v<remove_cvref_t<V>, InputType> &&
                       std::is_constructible_v<remove_cvref_t<V>, OutputType>,
                   remove_cvref_t<V>> {
    using EvalType = remove_cvref_t<V>;
    const auto xi = mapFromDomain(EvalType(pt));
    return detail::hybrid_impl<NCOEFFS, EvalType>(xi, coeffsBuf.data(), coeffsBuf.size());
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<int OuterUnrollFactor, bool ptsAligned, bool outAligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, NCOEFFS, ITERS, FUSION>::evalBatch(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
    return horner<NCOEFFS, ptsAligned, outAligned, OuterUnrollFactor>(
        pts, out, numPoints, coeffsBuf.data(), coeffsBuf.size(),
        [this](const auto v) { return this->mapFromDomain(v); });
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool ptsAligned, bool outAligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
#ifdef PF_OUTER_UNROLL
    PF_STATIC_CONSTEXPR_LOCAL auto unrollFactor = PF_OUTER_UNROLL;
#else
    PF_STATIC_CONSTEXPR_LOCAL auto unrollFactor = 0;
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

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR const typename FuncEval<Func, NCOEFFS, ITERS, FUSION>::OutputBuffer &
FuncEval<Func, NCOEFFS, ITERS, FUSION>::coeffs() const noexcept {
    return coeffsBuf;
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
constexpr std::size_t FuncEval<Func, NCOEFFS, ITERS, FUSION>::nCoeffs() const noexcept {
    return coeffsBuf.size();
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<class T>
[[nodiscard]] PF_ALWAYS_INLINE constexpr T FuncEval<Func, NCOEFFS, ITERS, FUSION>::mapToDomain(
    const DomainParams &dp, const T value) noexcept {
    if (dp.identityMap) return value;
    return polyfit::internal::helpers::mapToDomainScalar(value, dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<class T>
[[nodiscard]] PF_ALWAYS_INLINE constexpr T FuncEval<Func, NCOEFFS, ITERS, FUSION>::mapFromDomain(
    const T value) const noexcept {
    if constexpr (FUSION == FusionMode::Always) {
        return value;
    } else {
        if (domain_.identityMap) return value;
        return polyfit::internal::helpers::mapFromDomainScalar(value, domain_.invSpan, domain_.sumEndpoints);
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::initializeCoeffs(Func F, const InputType *pts,
                                                                                DomainParams &dp) {
    auto grid = makeBuffer<InputType, NCOEFFS>(coeffsBuf.size());
    auto samples = makeBuffer<OutputType, NCOEFFS>(coeffsBuf.size());

    buildNodeGrid(grid, pts);
    sampleOnGrid(samples, grid, F, dp);
    computeMonomialCoeffs(grid, samples);
    refine(grid, samples);

    if constexpr (FUSION != FusionMode::Never) {
        if (shouldFuseDomain(dp)) fuseDomain(dp);
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void
FuncEval<Func, NCOEFFS, ITERS, FUSION>::buildNodeGrid(InputBuffer &grid, const InputType *pts) const {
    const auto coeffCount = coeffsBuf.size();
    for (std::size_t coeffIdx = 0; coeffIdx < coeffCount; ++coeffIdx) {
        grid[coeffIdx] = pts ? pts[coeffIdx]
                             : InputType(detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi /
                                                     (2.0 * double(coeffCount))));
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::sampleOnGrid(
    OutputBuffer &samples, const InputBuffer &grid, Func F, const DomainParams &dp) const {
    const auto coeffCount = coeffsBuf.size();
    for (std::size_t sampleIdx = 0; sampleIdx < coeffCount; ++sampleIdx) {
        samples[sampleIdx] = F(mapToDomain(dp, grid[sampleIdx]));
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::computeMonomialCoeffs(
    const InputBuffer &grid, const OutputBuffer &samples) {
    auto newtonCoeffs = detail::bjorckPereyra<NCOEFFS, InputType, OutputType>(grid, samples);
    auto monomialCoeffs = detail::newtonToMonomial<NCOEFFS, InputType, OutputType>(newtonCoeffs, grid);
    assert(monomialCoeffs.size() == coeffsBuf.size() && "size mismatch!");
    std::copy(monomialCoeffs.begin(), monomialCoeffs.end(), coeffsBuf.begin());
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR bool FuncEval<Func, NCOEFFS, ITERS, FUSION>::shouldFuseDomain(const DomainParams &dp) const noexcept {
    if constexpr (FUSION == FusionMode::Always) {
        return true;
    } else {
        using Scalar = detail::value_type_or_t<InputType>;
        const auto alpha = Scalar(2) * static_cast<Scalar>(dp.invSpan);
        const auto beta = -static_cast<Scalar>(dp.sumEndpoints) * static_cast<Scalar>(dp.invSpan);
        const auto coeffCount = static_cast<int>(coeffsBuf.size());
        const auto condBase = detail::math::abs(alpha) + detail::math::abs(beta) + Scalar(1);
        constexpr auto maxLog = Scalar(std::numeric_limits<Scalar>::digits10 - 3);
        return coeffCount > 1 && Scalar(coeffCount - 1) * detail::math::log10(condBase) < maxLog;
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::fuseDomain(DomainParams &dp) {
    using Scalar = detail::value_type_or_t<InputType>;
    const auto alpha = Scalar(2) * static_cast<Scalar>(dp.invSpan);
    const auto beta = -static_cast<Scalar>(dp.sumEndpoints) * static_cast<Scalar>(dp.invSpan);

    polyfit::internal::helpers::fuseLinearMap(coeffsBuf.data(), coeffsBuf.size(), alpha, beta);
    dp.invSpan = InputType(0.5);
    dp.sumEndpoints = InputType(0);
    dp.identityMap = true;
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void
FuncEval<Func, NCOEFFS, ITERS, FUSION>::refine(const InputBuffer &chebNodes, const OutputBuffer &samples) {
    const auto nCoeffs = coeffsBuf.size();
    std::reverse(coeffsBuf.begin(), coeffsBuf.end());

    const std::size_t totalIters = ITERS + (nCoeffs > 32 ? 2 : 0);
    for (std::size_t pass = 0; pass < totalIters; ++pass) {
        auto residuals = makeBuffer<OutputType, NCOEFFS>(nCoeffs);
        for (std::size_t sampleIdx = 0; sampleIdx < nCoeffs; ++sampleIdx) {
            auto polyValue = poly_eval::compensated_horner<NCOEFFS>(chebNodes[sampleIdx], coeffsBuf.data(), nCoeffs);
            residuals[sampleIdx] = samples[sampleIdx] - polyValue;
        }
        auto newtonResidual = detail::bjorckPereyra<NCOEFFS, InputType, OutputType>(chebNodes, residuals);
        auto monomialResidual = detail::newtonToMonomial<NCOEFFS, InputType, OutputType>(newtonResidual, chebNodes);

        for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx) {
            coeffsBuf[nCoeffs - 1 - coeffIdx] += monomialResidual[coeffIdx];
        }
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_CXX20_CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::truncate(
    detail::value_type_or_t<OutputType> eps) noexcept {
    if constexpr (NCOEFFS == 0) {
        std::size_t skip = 0;
        while (skip + 1 < coeffsBuf.size() && std::abs(coeffsBuf[skip]) < eps) ++skip;
        if (skip > 0) coeffsBuf.erase(coeffsBuf.begin(), coeffsBuf.begin() + static_cast<std::ptrdiff_t>(skip));
    }
}

} // namespace poly_eval

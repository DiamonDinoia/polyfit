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

constexpr int validatePositiveCoeffCount(const int n) {
    if (n <= 0) {
        throw std::invalid_argument("nCoeffs must be positive");
    }
    return n;
}

template<typename T> constexpr void validateDomain(const T &a, const T &b) {
    if constexpr (detail::hasTupleSize_v<T>) {
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

template<class Evaluator, class Func, class Points>
[[nodiscard]] double maxRelativeError(Func F, const Evaluator &evaluator, const Points &points) {
    double maxErr = 0.0;
    for (const auto &pt : points) {
        maxErr = std::max(maxErr, relativeL2Norm(F(pt), evaluator(pt)));
    }
    return maxErr;
}

template<class Evaluator, class Func, class Spec, class InputType>
[[nodiscard]] Evaluator fitToTolerance(Func F, Spec tolerance, InputType a, InputType b, std::size_t evalPointCount,
                                       std::size_t maxCoeffCount) {
    if (tolerance <= Spec(0)) {
        throw std::invalid_argument("Requested error tolerance must be positive");
    }

    validateDomain(a, b);
    const auto evalPoints = linspace(a, b, int(evalPointCount));

    for (std::size_t candidateCoeffCount = 1; candidateCoeffCount <= maxCoeffCount; ++candidateCoeffCount) {
        Evaluator evaluator(F, int(candidateCoeffCount), a, b);
        if (maxRelativeError(F, evaluator, evalPoints) <= tolerance) {
            return evaluator;
        }
    }

    const Evaluator fallback(F, int(maxCoeffCount), a, b);
    throw std::runtime_error("No coefficient count found for requested error tolerance. eps=" +
                             std::to_string(tolerance) + ", maxCoeffCount=" + std::to_string(maxCoeffCount) +
                             ", maxErr=" + std::to_string(maxRelativeError(F, fallback, evalPoints)));
}

} // namespace detail

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR FuncEval<Func, NCOEFFS, ITERS, FUSION>::FuncEval(Func F, const int nCoeffs, const InputType a,
                                                                 const InputType b, const InputType *pts)
    : invSpan(InputType(1) / (b - a)), sumEndpoints(b + a) {
    initialize(detail::RuntimeCountTag{}, F, nCoeffs, a, b, pts);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR FuncEval<Func, NCOEFFS, ITERS, FUSION>::FuncEval(Func F, const InputType a, const InputType b,
                                                                 const InputType *pts)
    : invSpan(InputType(1) / (b - a)), sumEndpoints(b + a) {
    initialize(detail::CompileTimeCountTag{}, F, a, b, pts);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void
FuncEval<Func, NCOEFFS, ITERS, FUSION>::initialize(detail::CompileTimeCountTag, Func F, const InputType a,
                                                   const InputType b, const InputType *pts) {
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    detail::validateDomain(a, b);
    initializeCoeffs(F, pts);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::initialize(detail::RuntimeCountTag, Func F,
                                                                         const int nCoeffs, const InputType a,
                                                                         const InputType b, const InputType *pts) {
    static_assert(NCOEFFS == 0, "Runtime coefficient count is only valid for runtime-sized evaluators");
    detail::validateDomain(a, b);
    const auto validatedNCoeffs = detail::validatePositiveCoeffCount(nCoeffs);
    coeffsBuf.resize(static_cast<std::size_t>(validatedNCoeffs));
    initializeCoeffs(F, pts);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool>
constexpr typename FuncEval<Func, NCOEFFS, ITERS, FUSION>::OutputType PF_ALWAYS_INLINE
FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(const InputType pt) const noexcept {
    const auto xi = mapFromDomain(pt);
    return horner<NCOEFFS>(xi, coeffsBuf.data(), coeffsBuf.size());
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

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR const typename FuncEval<Func, NCOEFFS, ITERS, FUSION>::OutputBuffer &
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
    const T value) const noexcept {
    return polyfit::internal::helpers::mapToDomainScalar(value, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<class T>
[[nodiscard]] PF_ALWAYS_INLINE constexpr T FuncEval<Func, NCOEFFS, ITERS, FUSION>::mapFromDomain(
    const T value) const noexcept {
    if constexpr (FUSION == FusionMode::Always)
        return value;
    else if constexpr (FUSION == FusionMode::Never)
        return polyfit::internal::helpers::mapFromDomainScalar(value, invSpan, sumEndpoints);
    else {
        if (domainFused) return value;
        return polyfit::internal::helpers::mapFromDomainScalar(value, invSpan, sumEndpoints);
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::initializeCoeffs(Func F, const InputType *pts) {
    auto grid = makeBuffer<InputType, NCOEFFS>(coeffsBuf.size());
    auto samples = makeBuffer<OutputType, NCOEFFS>(coeffsBuf.size());

    buildNodeGrid(grid, pts);
    sampleOnGrid(samples, grid, F);
    computeMonomialCoeffs(grid, samples);
    refine(grid, samples);

    if constexpr (FUSION != FusionMode::Never) {
#if __cplusplus >= 202602L
        {
#else
        PF_IF_NOT_CONSTEVAL {
#endif
            if (shouldFuseDomain()) fuseDomain();
        }
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void
FuncEval<Func, NCOEFFS, ITERS, FUSION>::buildNodeGrid(InputBuffer &grid, const InputType *pts) const {
    const auto coeffCount = coeffsBuf.size();
    for (std::size_t coeffIdx = 0; coeffIdx < coeffCount; ++coeffIdx) {
        grid[coeffIdx] = pts ? pts[coeffIdx]
                             : InputType(detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi /
                                                     (2.0 * double(coeffCount))));
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::sampleOnGrid(
    OutputBuffer &samples, const InputBuffer &grid, Func F) const {
    const auto coeffCount = coeffsBuf.size();
    for (std::size_t sampleIdx = 0; sampleIdx < coeffCount; ++sampleIdx) {
        samples[sampleIdx] = F(mapToDomain(grid[sampleIdx]));
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::computeMonomialCoeffs(
    const InputBuffer &grid, const OutputBuffer &samples) {
    auto newtonCoeffs = detail::bjorckPereyra<NCOEFFS, InputType, OutputType>(grid, samples);
    auto monomialCoeffs = detail::newtonToMonomial<NCOEFFS, InputType, OutputType>(newtonCoeffs, grid);
    assert(monomialCoeffs.size() == coeffsBuf.size() && "size mismatch!");
    std::copy(monomialCoeffs.begin(), monomialCoeffs.end(), coeffsBuf.begin());
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR bool FuncEval<Func, NCOEFFS, ITERS, FUSION>::shouldFuseDomain() const noexcept {
    if constexpr (FUSION == FusionMode::Always) {
        return true;
    } else {
        using Scalar = typename detail::valueTypeOrIdentity<InputType>::type;
        const auto alpha = Scalar(2) * static_cast<Scalar>(invSpan);
        const auto beta = -static_cast<Scalar>(sumEndpoints) * static_cast<Scalar>(invSpan);
        const auto coeffCount = static_cast<int>(coeffsBuf.size());
        const auto condBase = detail::math::abs(alpha) + detail::math::abs(beta) + Scalar(1);
        constexpr auto maxLog = Scalar(std::numeric_limits<Scalar>::digits10 - 3);
        return coeffCount > 1 && Scalar(coeffCount - 1) * detail::math::log10(condBase) < maxLog;
    }
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::fuseDomain() {
    using Scalar = typename detail::valueTypeOrIdentity<InputType>::type;
    const auto alpha = Scalar(2) * static_cast<Scalar>(invSpan);
    const auto beta = -static_cast<Scalar>(sumEndpoints) * static_cast<Scalar>(invSpan);

    polyfit::internal::helpers::fuseLinearMap(coeffsBuf.data(), coeffsBuf.size(), alpha, beta);
    invSpan = InputType(0.5);
    sumEndpoints = InputType(0);
    if constexpr (FUSION == FusionMode::Auto) domainFused = true;
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
PF_C20CONSTEXPR void
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
PF_C20CONSTEXPR void FuncEval<Func, NCOEFFS, ITERS, FUSION>::truncate(
    typename detail::valueTypeOrIdentity<OutputType>::type eps) noexcept {
    if constexpr (NCOEFFS == 0) {
        std::size_t skip = 0;
        while (skip + 1 < coeffsBuf.size() && std::abs(coeffsBuf[skip]) < eps) ++skip;
        if (skip > 0) coeffsBuf.erase(coeffsBuf.begin(), coeffsBuf.begin() + static_cast<std::ptrdiff_t>(skip));
    }
}

template<typename... EvalTypes> PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const EvalTypes &...evals) {
    invSpan = {evals.invSpan...};
    sumEndpoints = {evals.sumEndpoints...};

    if constexpr (MAX_NCOEFFS == 0) {
        bindCoeffView(std::max({std::size_t(evals.coeffsBuf.size())...}));
    } else {
        bindCoeffView(MAX_NCOEFFS);
    }

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

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::size() const noexcept { return COUNT; }

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::nCoeffs() const noexcept {
    return coeffs.extent(0);
}

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(InputType x) const noexcept -> std::array<OutputType, COUNT> {
    return evalMapped(mapInputs(x));
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(const std::array<InputType, COUNT> &xs) const noexcept
    -> std::array<OutputType, COUNT> {
    return evalMapped(mapInputs(xs));
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::operator()(const InputType *PF_RESTRICT x, OutputType *PF_RESTRICT out,
                                            std::size_t numPoints) const noexcept {
    poet::static_for<COUNT>([&](auto column) { evalColumn(std::size_t(column), x, out, numPoints); });
}
PF_FAST_EVAL_END

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(InputType first, Ts... rest) const noexcept -> std::array<OutputType, COUNT> {
    static_assert(sizeof...(Ts) + 1 == COUNT, "Incorrect number of arguments");
    return operator()(std::array<InputType, COUNT>{first, static_cast<InputType>(rest)...});
}

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(const std::tuple<Ts...> &tup) const noexcept -> std::array<OutputType, COUNT> {
    static_assert(sizeof...(Ts) == COUNT, "Tuple size must equal number of polynomials");
    std::array<InputType, COUNT> xs{};
    std::apply([&](auto &&...e) { xs = {static_cast<InputType>(e)...}; }, tup);
    return operator()(xs);
}

template<typename... EvalTypes>
template<std::size_t I, typename FE, typename... Rest>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::copyCoeffs(const FE &eval, const Rest &...rest) {
    for (std::size_t k = 0; k < eval.coeffsBuf.size(); ++k) coeffs(k, I) = eval.coeffsBuf[k];
    for (std::size_t k = eval.coeffsBuf.size(); k < coeffs.extent(0); ++k) coeffs(k, I) = OutputType{};
    if constexpr (I + 1 < COUNT) copyCoeffs<I + 1>(rest...);
}

template<typename... EvalTypes>
template<class Step>
constexpr void FuncEvalMany<EvalTypes...>::forEachCoeff(Step &&step) const noexcept {
    if constexpr (MAX_NCOEFFS != 0) {
        poet::static_for<1, MAX_NCOEFFS>(std::forward<Step>(step));
    } else {
        for (std::size_t k = 1; k < coeffs.extent(0); ++k) step(k);
    }
}

template<typename... EvalTypes> PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::bindCoeffView(std::size_t coeffCount) {
    if constexpr (MAX_NCOEFFS == 0) {
        coeffStore.assign(PADDED_COUNT * coeffCount, OutputType{});
    }
    coeffs = decltype(coeffs){coeffStore.data(), coeffCount, PADDED_COUNT};
}

template<typename... EvalTypes>
constexpr typename FuncEvalMany<EvalTypes...>::InputType
FuncEvalMany<EvalTypes...>::mapInput(std::size_t polyIndex, InputType x) const noexcept {
    return xsimd::fms(InputType(2.0), x, sumEndpoints[polyIndex]) * invSpan[polyIndex];
}

template<typename... EvalTypes>
constexpr std::array<typename FuncEvalMany<EvalTypes...>::InputType, FuncEvalMany<EvalTypes...>::PADDED_COUNT>
FuncEvalMany<EvalTypes...>::mapInputs(InputType x) const noexcept {
    std::array<InputType, PADDED_COUNT> mapped{};
    poet::static_for<COUNT>([&](auto i) { mapped[i] = mapInput(std::size_t(i), x); });
    return mapped;
}

template<typename... EvalTypes>
constexpr std::array<typename FuncEvalMany<EvalTypes...>::InputType, FuncEvalMany<EvalTypes...>::PADDED_COUNT>
FuncEvalMany<EvalTypes...>::mapInputs(const std::array<InputType, COUNT> &xs) const noexcept {
    std::array<InputType, PADDED_COUNT> mapped{};
    poet::static_for<COUNT>([&](auto i) { mapped[i] = mapInput(std::size_t(i), xs[i]); });
    return mapped;
}

template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::evalMapped(const std::array<InputType, PADDED_COUNT> &xu) const noexcept
    -> std::array<OutputType, COUNT> {
    alignas(ALIGNMENT) std::array<OutputType, PADDED_COUNT> full{};
    horner_transposed<PADDED_COUNT, MAX_NCOEFFS, VECTOR_WIDTH, true>(
        xu.data(), coeffs.data_handle(), full.data(), PADDED_COUNT, static_cast<std::size_t>(coeffs.extent(0)));
    return extractReal(full);
}

template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::scatterColumnBatch(xsimd::batch<OutputType> acc, OutputType *out, std::size_t base,
                                                    std::size_t column) const noexcept {
    constexpr std::size_t simdSize = xsimd::batch<InputType>::size;
    alignas(xsimd::batch<OutputType>::arch_type::alignment()) OutputType tmp[simdSize];
    acc.store_aligned(tmp);
    poet::static_for<simdSize>([&](auto lane) { out[(base + std::size_t(lane)) * COUNT + column] = tmp[lane]; });
}

template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::evalColumn(std::size_t column, const InputType *PF_RESTRICT x,
                                            OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
    PF_C23STATIC constexpr std::size_t simdSize = xsimd::batch<InputType>::size;
    PF_C23STATIC constexpr std::size_t unrollFactor = detail::optimalManyEvalUf<OutputType>();
    PF_C23STATIC constexpr std::size_t blockSize = simdSize * unrollFactor;
    PF_C23STATIC constexpr std::size_t stride = PADDED_COUNT;

    const auto invSpanValue = invSpan[column];
    const auto sumEndpointsValue = sumEndpoints[column];
    const OutputType *columnCoeffs = coeffs.data_handle() + column;
    const auto twoVec = xsimd::batch<InputType>(InputType(2.0));
    const auto sumEndpointsVec = xsimd::batch<InputType>(sumEndpointsValue);
    const auto invSpanVec = xsimd::batch<InputType>(invSpanValue);
    auto mapSimd = [&](auto xv) { return xsimd::fms(twoVec, xv, sumEndpointsVec) * invSpanVec; };

    const auto tileEnd = detail::roundDown<blockSize>(numPoints);
    const auto simdEnd = detail::roundDown<simdSize>(numPoints);

    for (std::size_t i = 0; i < tileEnd; i += blockSize) {
        xsimd::batch<InputType> pt[unrollFactor];
        xsimd::batch<OutputType> acc[unrollFactor];

        poet::static_for<unrollFactor>([&](auto j) {
            pt[j] = mapSimd(xsimd::load_unaligned(x + i + j * simdSize));
            acc[j] = xsimd::batch<OutputType>(columnCoeffs[0]);
        });

        forEachCoeff([&](auto k) {
            const auto coeff = xsimd::batch<OutputType>(columnCoeffs[std::size_t(k) * stride]);
            poet::static_for<unrollFactor>([&](auto j) { acc[j] = detail::fma(acc[j], pt[j], coeff); });
        });

        poet::static_for<unrollFactor>([&](auto j) { scatterColumnBatch(acc[j], out, i + j * simdSize, column); });
    }

    for (std::size_t i = tileEnd; i < simdEnd; i += simdSize) {
        auto mapped = mapSimd(xsimd::load_unaligned(x + i));
        auto acc = xsimd::batch<OutputType>(columnCoeffs[0]);
        forEachCoeff([&](auto k) {
            acc = detail::fma(acc, mapped, xsimd::batch<OutputType>(columnCoeffs[std::size_t(k) * stride]));
        });
        scatterColumnBatch(acc, out, i, column);
    }

    for (std::size_t i = simdEnd; i < numPoints; ++i) {
        auto mapped = mapInput(column, x[i]);
        OutputType acc = columnCoeffs[0];
        forEachCoeff([&](auto k) { acc = detail::fma(acc, mapped, columnCoeffs[std::size_t(k) * stride]); });
        out[i * COUNT + column] = acc;
    }
}

template<typename... EvalTypes> PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::zeroPadCoeffs() {
    for (std::size_t j = COUNT; j < PADDED_COUNT; ++j)
        for (std::size_t k = 0; k < coeffs.extent(0); ++k) coeffs(k, j) = OutputType{};
}

template<typename... EvalTypes>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::truncate(typename detail::valueTypeOrIdentity<OutputType>::type eps) {
    std::size_t newNCoeffs = coeffs.extent(0);
    while (newNCoeffs > 1) {
        OutputType rowMax{};
        for (std::size_t j = 0; j < COUNT; ++j) rowMax = std::max(rowMax, std::abs(coeffs(newNCoeffs - 1, j)));
        if (rowMax >= eps) break;
        --newNCoeffs;
    }
    if (newNCoeffs < coeffs.extent(0)) {
        coeffs = decltype(coeffs){coeffStore.data(), newNCoeffs, PADDED_COUNT};
    }
}

template<typename... EvalTypes>
constexpr auto FuncEvalMany<EvalTypes...>::extractReal(const std::array<OutputType, PADDED_COUNT> &full) const noexcept
    -> std::array<OutputType, COUNT> {
    if constexpr (COUNT == PADDED_COUNT) {
        return full;
    }
    std::array<OutputType, COUNT> out{};
    poet::static_for<COUNT>([&](auto i) { out[i] = full[i]; });
    return out;
}

template<class Func, std::size_t NCOEFFS>
constexpr FuncEvalND<Func, NCOEFFS>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : coeffsFlat(), coeffsMd{coeffsFlat.data(), Extents{}} {
    initialize(detail::CompileTimeCountTag{}, f, a, b);
}

template<class Func, std::size_t NCOEFFS>
constexpr FuncEvalND<Func, NCOEFFS>::FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b)
    : coeffsFlat(storageRequired(detail::validatePositiveCoeffCount(nCoeffsPerAxis))),
      coeffsMd{coeffsFlat.data(), makeExtents(nCoeffsPerAxis)} {
    initialize(detail::RuntimeCountTag{}, f, nCoeffsPerAxis, a, b);
}

template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::initialize(detail::CompileTimeCountTag, Func f, const InputType &a,
                                                     const InputType &b) {
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    detail::validateDomain(a, b);
    computeScaling(a, b);
    buildCoeffs(static_cast<int>(NCOEFFS), f);
}

template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::initialize(detail::RuntimeCountTag, Func f, int nCoeffsPerAxis,
                                                     const InputType &a, const InputType &b) {
    static_assert(NCOEFFS == 0, "Runtime coefficient count is only valid for runtime-sized evaluators");
    detail::validateDomain(a, b);
    computeScaling(a, b);
    buildCoeffs(nCoeffsPerAxis, f);
}

template<class Func, std::size_t NCOEFFS>
constexpr std::size_t FuncEvalND<Func, NCOEFFS>::nCoeffsPerAxis() const noexcept {
    return static_cast<std::size_t>(coeffsMd.extent(0));
}

template<class Func, std::size_t NCOEFFS>
FuncEvalND<Func, NCOEFFS>::FuncEvalND(const FuncEvalND &other)
    : invSpan(other.invSpan), sumEndpoints(other.sumEndpoints), coeffsFlat(other.coeffsFlat),
      coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS>
auto FuncEvalND<Func, NCOEFFS>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        invSpan = other.invSpan;
        sumEndpoints = other.sumEndpoints;
        coeffsFlat = other.coeffsFlat;
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

template<class Func, std::size_t NCOEFFS>
FuncEvalND<Func, NCOEFFS>::FuncEvalND(FuncEvalND &&other) noexcept
    : invSpan(std::move(other.invSpan)), sumEndpoints(std::move(other.sumEndpoints)),
      coeffsFlat(std::move(other.coeffsFlat)), coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS>
auto FuncEvalND<Func, NCOEFFS>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        invSpan = std::move(other.invSpan);
        sumEndpoints = std::move(other.sumEndpoints);
        coeffsFlat = std::move(other.coeffsFlat);
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCOEFFS>::OutputType FuncEvalND<Func, NCOEFFS>::operator()(
    const InputType &x) const {
    const int nCoeffsRt = (NCOEFFS ? static_cast<int>(NCOEFFS) : static_cast<int>(coeffsMd.extent(0)));
    return poly_eval::horner<NCOEFFS, SIMD, OutputType>(mapFromDomain(x), coeffsMd, nCoeffsRt);
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCOEFFS>
template<typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, NCOEFFS>::Scalar &FuncEvalND<Func, NCOEFFS>::coeffImpl(
    const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffsMd(static_cast<std::size_t>(idx[I])..., k);
}

template<class Func, std::size_t NCOEFFS>
template<class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS>::Scalar &FuncEvalND<Func, NCOEFFS>::coeff(
    const IdxArray &idx, std::size_t k) noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<DIM>{});
}

template<class Func, std::size_t NCOEFFS>
auto FuncEvalND<Func, NCOEFFS>::makeExtents(int nCoeffsPerAxis) noexcept -> Extents {
    if constexpr (IS_STATIC) {
        return detail::makeStaticExtents<NCOEFFS, DIM, OUT_DIM>(std::make_index_sequence<DIM>{});
    } else {
        return makeExtents(nCoeffsPerAxis, std::make_index_sequence<DIM + 1>{});
    }
}

template<class Func, std::size_t NCOEFFS>
template<std::size_t... Is>
auto FuncEvalND<Func, NCOEFFS>::makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept -> Extents {
    return Extents{(Is < DIM ? static_cast<std::size_t>(nCoeffsPerAxis) : static_cast<std::size_t>(OUT_DIM))...};
}

template<class Func, std::size_t NCOEFFS>
constexpr std::size_t FuncEvalND<Func, NCOEFFS>::storageRequired(const int nCoeffsPerAxis) noexcept {
    auto ext = makeExtents(nCoeffsPerAxis);
    auto mapping = typename Mdspan::mapping_type{ext};
    return mapping.required_span_size();
}
template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::buildCoeffs(int nCoeffsPerAxis, Func f) {
    const auto nCoeffs = static_cast<std::size_t>(nCoeffsPerAxis);
    auto nodes = makeBuffer<Scalar, NCOEFFS>(nCoeffs);
    for (std::size_t coeffIdx = 0; coeffIdx < nCoeffs; ++coeffIdx)
        nodes[coeffIdx] =
            detail::cos((2.0 * double(coeffIdx) + 1.0) * detail::constants::pi / (2.0 * double(nCoeffsPerAxis)));

    std::array<int, DIM> extents{};
    extents.fill(nCoeffsPerAxis);

    // sample f on Chebyshev grid
    forEachIndex<DIM>(extents, [&](const std::array<int, DIM> &idx) {
        InputType domainPoint{};
        for (std::size_t d = 0; d < DIM; ++d) domainPoint[d] = nodes[static_cast<std::size_t>(idx[d])];
        OutputType y = f(mapToDomain(domainPoint));
        for (std::size_t k = 0; k < OUT_DIM; ++k) coeff(idx, k) = y[k];
    });

    convertAxesToMonomialBasis(nCoeffsPerAxis, nodes);
    reverseCoefficientOrder(nCoeffsPerAxis);
}

template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::convertAxesToMonomialBasis(int nCoeffsPerAxis,
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

template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::reverseCoefficientOrder(int nCoeffsPerAxis) {
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

template<class Func, std::size_t NCOEFFS>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS>::InputType FuncEvalND<Func, NCOEFFS>::mapToDomain(
    const InputType &x) const noexcept {
    return polyfit::internal::helpers::mapToDomainArray<Scalar, DIM>(x, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCOEFFS>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS>::InputType FuncEvalND<Func, NCOEFFS>::mapFromDomain(
    const InputType &x) const noexcept {
    return polyfit::internal::helpers::mapFromDomainArray<Scalar, DIM>(x, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCOEFFS>
constexpr void FuncEvalND<Func, NCOEFFS>::computeScaling(const InputType &a, const InputType &b) noexcept {
    polyfit::internal::helpers::computeScalingArray<Scalar, DIM>(a, b, invSpan, sumEndpoints);
}

template<class Func, std::size_t NCOEFFS>
template<std::size_t Rank, class F>
constexpr void FuncEvalND<Func, NCOEFFS>::forEachIndex(const std::array<int, Rank> &ext, F &&body) {
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

template<std::size_t NCOEFFS, class Func, class... Tags>
[[nodiscard]] PF_C20CONSTEXPR auto fit(Func F, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported fit tag");
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    using Evaluator = FitEvaluator<Func, NCOEFFS, Options::ITERS, Options::FUSION_MODE>;
    return Evaluator(F, a, b);
}

template<class Func, class Spec, class... Tags>
[[nodiscard]] PF_C20CONSTEXPR auto fit(Func F, Spec spec, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported fit tag");

    if constexpr (isIntegralLike_v<Spec>) {
        const auto nCoeffs = detail::validatePositiveCoeffCount(static_cast<int>(spec));
        if constexpr (takesTupleInput_v<Func>) {
            return FuncEvalND<Func, 0>(F, nCoeffs, a, b);
        } else {
            using Evaluator = FitEvaluator<Func, 0, Options::ITERS, Options::FUSION_MODE>;
            return Evaluator(F, nCoeffs, a, b);
        }
    } else if constexpr (isFloatingPointLike_v<Spec>) {
        using Evaluator = FitEvaluator<Func, 0, Options::ITERS, Options::FUSION_MODE>;
        return detail::fitToTolerance<Evaluator>(F, spec, a, b, Options::EVAL_POINTS, Options::MAX_NCOEFFS);
    } else {
        static_assert(alwaysFalse_v<Spec>,
                      "fit(...) expects an integral coefficient count or a floating-point tolerance");
    }
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double EPS, auto a, auto b, std::size_t MAX_NCOEFFS, std::size_t EVAL_POINTS, std::size_t ITERS, class Func>
[[nodiscard]] constexpr auto fit(Func F) {
    using RawInputType = poly_eval::remove_cvref_t<typename FunctionTraits<Func>::arg0_type>;
    static_assert(MAX_NCOEFFS > 0, "Max coefficient count must be positive.");
    static_assert(EVAL_POINTS > 1, "Number of evaluation points must be greater than 1.");

    constexpr auto nCoeffs = [F] {
        constexpr auto computeError = [F](const auto &evaluator) {
            constexpr auto ep = detail::linspace<static_cast<int>(EVAL_POINTS)>(a, b);
            double maxErr = 0.0;
            for (const auto &pt : ep) {
                const auto actual = F(pt);
                const auto approx = evaluator.template operator()<false>(pt);
                maxErr = std::max(detail::relativeL2Norm(actual, approx), maxErr);
            }
            return maxErr;
        };
        int result = 0;
        poet::static_for<1, MAX_NCOEFFS + 1>([&](auto i) {
            if (result != 0) return;
            using Evaluator = std::conditional_t<detail::hasTupleSize_v<RawInputType>, FuncEvalND<Func, i>,
                                                 FuncEval<Func, i, ITERS>>;
            if constexpr (computeError(Evaluator(F, a, b)) <= EPS) {
                result = i;
            }
        });
        return result;
    }();
    static_assert(nCoeffs != 0, "No coefficient count found for requested error tolerance.");
    using Evaluator = std::conditional_t<detail::hasTupleSize_v<RawInputType>, FuncEvalND<Func, nCoeffs>,
                                         FuncEval<Func, nCoeffs, ITERS>>;
    return Evaluator(F, a, b);
}
#endif

template<typename... EvalTypes>
[[nodiscard]] PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> pack(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval

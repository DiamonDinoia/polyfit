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
    if (n <= 0) PF_UNLIKELY {
        throw std::invalid_argument("nCoeffs must be positive");
    }
    return n;
}

template<typename T> constexpr void validateDomain(const T &a, const T &b) {
    if constexpr (detail::isFixedIndexable_v<T>) {
        using Access = detail::FixedContainerAccess<T>;
        for (std::size_t i = 0; i < Access::size; ++i) {
            if (Access::get(a, i) == Access::get(b, i)) PF_UNLIKELY {
                throw std::invalid_argument("Domain endpoints must differ in every dimension");
            }
        }
    } else {
        if (a == b) PF_UNLIKELY {
            throw std::invalid_argument("Domain endpoints must differ");
        }
    }
}

template<class Evaluator, class Func, class Points>
[[nodiscard]] double maxRelativeError(Func F, const Evaluator &evaluator, const Points &points) {
    double maxErr = 0.0;
    for (const auto &pt : points) {
        maxErr = std::max(maxErr, relativeL2Norm(evaluator(pt), F(pt)));
    }
    return maxErr;
}

template<class Evaluator, class Func, class Spec, class InputType>
[[nodiscard]] Evaluator fitToTolerance(Func F, Spec tolerance, InputType a, InputType b, std::size_t evalPointCount,
                                       std::size_t maxCoeffCount) {
    if (tolerance <= Spec(0)) PF_UNLIKELY {
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

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool>
constexpr typename FuncEval<Func, NCOEFFS, ITERS, FUSION>::OutputType PF_ALWAYS_INLINE
FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(const InputType pt) const noexcept {
    const auto xi = mapFromDomain(pt);
    return horner<NCOEFFS>(xi, coeffsBuf.data(), coeffsBuf.size());
}

template<class Func, std::size_t NCOEFFS, std::size_t ITERS, FusionMode FUSION>
template<bool Dummy, class V>
constexpr auto PF_ALWAYS_INLINE FuncEval<Func, NCOEFFS, ITERS, FUSION>::operator()(V pt) const noexcept
    -> enable_if_t<!std::is_same_v<remove_cvref_t<V>, InputType> &&
                       std::is_constructible_v<remove_cvref_t<V>, OutputType>,
                   remove_cvref_t<V>> {
    using EvalType = remove_cvref_t<V>;
    const auto xi = mapFromDomain(EvalType(pt));
    return detail::horner_impl<NCOEFFS, EvalType>(xi, coeffsBuf.data(), coeffsBuf.size());
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

template<typename... EvalTypes> PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const EvalTypes &...evals) {
    invSpan = {evals.domainInvSpan()...};
    sumEndpoints = {evals.domainSumEndpoints()...};
    identityMap = {evals.domainIsIdentity()...};
    allIdentityMap = (... && evals.domainIsIdentity());

    if constexpr (MAX_NCOEFFS == 0) {
        bindCoeffView(std::max({std::size_t(evals.coeffsBuf.size())...}));
    } else {
        bindCoeffView(MAX_NCOEFFS);
    }

    copyCoeffs<0>(evals...);
    zeroPadCoeffs();
}

template<typename... EvalTypes>
PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const FuncEvalMany &other)
    : coeffStore(other.coeffStore), coeffs{coeffStore.data(), other.coeffs.extents()}, invSpan(other.invSpan),
      sumEndpoints(other.sumEndpoints), identityMap(other.identityMap), allIdentityMap(other.allIdentityMap) {}

template<typename... EvalTypes>
PF_CXX20_CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(const FuncEvalMany &other) -> FuncEvalMany & {
    if (this != &other) {
        coeffStore = other.coeffStore;
        invSpan = other.invSpan;
        sumEndpoints = other.sumEndpoints;
        identityMap = other.identityMap;
        allIdentityMap = other.allIdentityMap;
        rebindCoeffs(other.coeffs.extent(0));
    }
    return *this;
}

template<typename... EvalTypes>
PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(FuncEvalMany &&other) noexcept
    : coeffStore(std::move(other.coeffStore)), coeffs{coeffStore.data(), other.coeffs.extents()},
      invSpan(std::move(other.invSpan)), sumEndpoints(std::move(other.sumEndpoints)), identityMap(other.identityMap),
      allIdentityMap(other.allIdentityMap) {}

template<typename... EvalTypes>
PF_CXX20_CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(FuncEvalMany &&other) noexcept -> FuncEvalMany & {
    if (this != &other) {
        coeffStore = std::move(other.coeffStore);
        invSpan = std::move(other.invSpan);
        sumEndpoints = std::move(other.sumEndpoints);
        identityMap = other.identityMap;
        allIdentityMap = other.allIdentityMap;
        rebindCoeffs(other.coeffs.extent(0));
    }
    return *this;
}

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::size() const noexcept { return COUNT; }

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::nCoeffs() const noexcept {
    return coeffs.extent(0);
}

template<typename... EvalTypes>
constexpr const typename FuncEvalMany<EvalTypes...>::OutputType &
FuncEvalMany<EvalTypes...>::coeff(std::size_t coeffIndex, std::size_t polyIndex) const noexcept {
    return coeffRef(coeffIndex, polyIndex);
}

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(InputType x) const noexcept -> std::array<OutputType, COUNT> {
    if (allIdentityMap) {
        std::array<InputType, PADDED_COUNT> inputs{};
        inputs.fill(x);
        return evalInputs(inputs);
    }
    return evalInputs(mapInputs(x));
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(const std::array<InputType, COUNT> &xs) const noexcept
    -> std::array<OutputType, COUNT> {
    if (allIdentityMap) return evalInputs(gatherInputs([&](std::size_t i) { return xs[i]; }));
    return evalInputs(mapInputs(xs));
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
PF_CXX20_CONSTEXPR void FuncEvalMany<EvalTypes...>::copyCoeffs(const FE &eval, const Rest &...rest) {
    for (std::size_t k = 0; k < eval.coeffsBuf.size(); ++k) coeffRef(k, I) = eval.coeffsBuf[k];
    for (std::size_t k = eval.coeffsBuf.size(); k < coeffs.extent(0); ++k) coeffRef(k, I) = OutputType{};
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

template<typename... EvalTypes> PF_CXX20_CONSTEXPR void FuncEvalMany<EvalTypes...>::rebindCoeffs(std::size_t coeffCount) noexcept {
    coeffs = decltype(coeffs){coeffStore.data(), coeffCount, PADDED_COUNT};
}

template<typename... EvalTypes> PF_CXX20_CONSTEXPR void FuncEvalMany<EvalTypes...>::bindCoeffView(std::size_t coeffCount) {
    if constexpr (MAX_NCOEFFS == 0) {
        coeffStore.assign(PADDED_COUNT * coeffCount, OutputType{});
    }
    rebindCoeffs(coeffCount);
}

template<typename... EvalTypes>
constexpr auto FuncEvalMany<EvalTypes...>::coeffRef(std::size_t coeffIndex, std::size_t polyIndex) noexcept -> OutputType & {
    return coeffs[std::array<std::size_t, 2>{coeffIndex, polyIndex}];
}

template<typename... EvalTypes>
constexpr auto FuncEvalMany<EvalTypes...>::coeffRef(std::size_t coeffIndex, std::size_t polyIndex) const noexcept
    -> const OutputType & {
    return coeffs[std::array<std::size_t, 2>{coeffIndex, polyIndex}];
}

template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::evalInputs(const std::array<InputType, PADDED_COUNT> &xu) const noexcept
    -> std::array<OutputType, COUNT> {
    alignas(ALIGNMENT) std::array<InputType, PADDED_COUNT> alignedXu = xu;
    alignas(ALIGNMENT) std::array<OutputType, PADDED_COUNT> full{};
    horner_transposed<PADDED_COUNT, MAX_NCOEFFS, VECTOR_WIDTH, true>(
        alignedXu.data(), coeffs.data_handle(), full.data(), PADDED_COUNT, static_cast<std::size_t>(coeffs.extent(0)));
    return extractReal(full);
}

template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::scatterColumnBatch(xsimd::batch<OutputType> acc, OutputType *out, std::size_t base,
                                                    std::size_t column) const noexcept {
    using Batch = xsimd::batch<OutputType>;
    constexpr std::size_t simdSize = Batch::size;
    alignas(Batch::arch_type::alignment()) OutputType tmp[simdSize];
    acc.store_aligned(tmp);
    poet::static_for<simdSize>([&](auto lane) { out[(base + std::size_t(lane)) * COUNT + column] = tmp[lane]; });
}

template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::evalColumn(std::size_t column, const InputType *PF_RESTRICT x,
                                            OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept {
    PF_STATIC_CONSTEXPR_LOCAL std::size_t stride = PADDED_COUNT;
    using InBatch = xsimd::batch<InputType>;
    using OutBatch = xsimd::batch<OutputType>;

    const auto invSpanValue = invSpan[column];
    const auto sumEndpointsValue = sumEndpoints[column];
    const OutputType *columnCoeffs = coeffs.data_handle() + column;
    const bool columnIdentity = identityMap[column];
    if constexpr (InBatch::size != OutBatch::size) {
        for (std::size_t i = 0; i < numPoints; ++i) {
            auto mapped = mapInput(column, x[i]);
            OutputType acc = columnCoeffs[0];
            forEachCoeff([&](auto k) { acc = detail::fma(acc, mapped, columnCoeffs[std::size_t(k) * stride]); });
            out[i * COUNT + column] = acc;
        }
        return;
    } else {
        PF_STATIC_CONSTEXPR_LOCAL std::size_t simdSize = OutBatch::size;
        PF_STATIC_CONSTEXPR_LOCAL std::size_t unrollFactor = detail::optimalManyEvalUf<OutputType>();
        PF_STATIC_CONSTEXPR_LOCAL std::size_t blockSize = simdSize * unrollFactor;
        const auto twoVec = InBatch(InputType(2.0));
        const auto sumEndpointsVec = InBatch(sumEndpointsValue);
        const auto invSpanVec = InBatch(invSpanValue);
        auto mapSimd = [&](InBatch xv) {
            if (columnIdentity) return xv;
            return xsimd::fms(twoVec, xv, sumEndpointsVec) * invSpanVec;
        };

        const auto tileEnd = detail::roundDown<blockSize>(numPoints);
        const auto simdEnd = detail::roundDown<simdSize>(numPoints);

        for (std::size_t i = 0; i < tileEnd; i += blockSize) {
            InBatch pt[unrollFactor];
            OutBatch acc[unrollFactor];

            poet::static_for<unrollFactor>([&](auto j) {
                pt[j] = mapSimd(xsimd::load_unaligned(x + i + j * simdSize));
                acc[j] = OutBatch(columnCoeffs[0]);
            });

            forEachCoeff([&](auto k) {
                const auto coeff = OutBatch(columnCoeffs[std::size_t(k) * stride]);
                poet::static_for<unrollFactor>([&](auto j) { acc[j] = detail::fma(acc[j], pt[j], coeff); });
            });

            poet::static_for<unrollFactor>([&](auto j) { scatterColumnBatch(acc[j], out, i + j * simdSize, column); });
        }

        for (std::size_t i = tileEnd; i < simdEnd; i += simdSize) {
            auto mapped = mapSimd(xsimd::load_unaligned(x + i));
            auto acc = OutBatch(columnCoeffs[0]);
            forEachCoeff([&](auto k) {
                acc = detail::fma(acc, mapped, OutBatch(columnCoeffs[std::size_t(k) * stride]));
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
}

template<typename... EvalTypes> PF_CXX20_CONSTEXPR void FuncEvalMany<EvalTypes...>::zeroPadCoeffs() {
    for (std::size_t j = COUNT; j < PADDED_COUNT; ++j)
        for (std::size_t k = 0; k < coeffs.extent(0); ++k) coeffRef(k, j) = OutputType{};
}

template<typename... EvalTypes>
PF_CXX20_CONSTEXPR void FuncEvalMany<EvalTypes...>::truncate(detail::value_type_or_t<OutputType> eps) {
    using Scalar = detail::value_type_or_t<OutputType>;
    std::size_t newNCoeffs = coeffs.extent(0);
    while (newNCoeffs > 1) {
        Scalar rowMax{};
        for (std::size_t j = 0; j < COUNT; ++j) {
            rowMax = std::max(rowMax, Scalar(std::abs(coeffRef(newNCoeffs - 1, j))));
        }
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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr FuncEvalND<Func, NCOEFFS, FUSION_MODE>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : coeffsFlat(), coeffsMd{coeffsFlat.data(), Extents{}} {
    static_assert(takesNdInput_v<Func>, "FuncEvalND requires fixed-size indexable ND input and output types");
    initialize(detail::CompileTimeCountTag{}, f, a, b);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr FuncEvalND<Func, NCOEFFS, FUSION_MODE>::FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b)
    : coeffsFlat(storageRequired(detail::validatePositiveCoeffCount(nCoeffsPerAxis))),
      coeffsMd{coeffsFlat.data(), makeExtents(nCoeffsPerAxis)} {
    static_assert(takesNdInput_v<Func>, "FuncEvalND requires fixed-size indexable ND input and output types");
    initialize(detail::RuntimeCountTag{}, f, nCoeffsPerAxis, a, b);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::initialize(detail::CompileTimeCountTag, Func f, const InputType &a,
                                                                  const InputType &b) {
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    detail::validateDomain(a, b);
    DomainParams dp;
    computeScaling(a, b, dp);
    buildCoeffs(static_cast<int>(NCOEFFS), f, dp);
    if constexpr (FUSION_MODE != FusionMode::Never) {
        fuseNDDomain(dp, static_cast<int>(NCOEFFS));
    }
    if constexpr (kStoresDomain) domain_ = dp;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::initialize(detail::RuntimeCountTag, Func f, int nCoeffsPerAxis,
                                                                  const InputType &a, const InputType &b) {
    static_assert(NCOEFFS == 0, "Runtime coefficient count is only valid for runtime-sized evaluators");
    detail::validateDomain(a, b);
    DomainParams dp;
    computeScaling(a, b, dp);
    buildCoeffs(nCoeffsPerAxis, f, dp);
    if constexpr (FUSION_MODE != FusionMode::Never) {
        fuseNDDomain(dp, nCoeffsPerAxis);
    }
    if constexpr (kStoresDomain) domain_ = dp;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr std::size_t FuncEvalND<Func, NCOEFFS, FUSION_MODE>::nCoeffsPerAxis() const noexcept {
    return static_cast<std::size_t>(coeffsMd.extent(0));
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::FuncEvalND(const FuncEvalND &other)
    : domain_(other.domain_), coeffsFlat(other.coeffsFlat),
      coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        domain_ = other.domain_;
        coeffsFlat = other.coeffsFlat;
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::FuncEvalND(FuncEvalND &&other) noexcept
    : domain_(std::move(other.domain_)),
      coeffsFlat(std::move(other.coeffsFlat)), coeffsMd{coeffsFlat.data(), other.coeffsMd.extents()} {}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        domain_ = std::move(other.domain_);
        coeffsFlat = std::move(other.coeffsFlat);
        coeffsMd = Mdspan{coeffsFlat.data(), other.coeffsMd.extents()};
    }
    return *this;
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(
    const InputType &x) const {
    return evalPoint<SIMD>(x);
}
PF_FAST_EVAL_END

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<bool SIMD, class Point, class>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(const Point &x) const {
    return evalPoint<SIMD>(x);
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<class... Coords, class>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::OutputType FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(Coords... coords) const {
    return fromCanonicalOutput(evalCanonical(CanonicalInput{static_cast<InputScalar>(coords)...}));
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<bool SIMD, class Point>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::OutputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::evalPoint(const Point &x) const {
    return fromCanonicalOutput(evalCanonical<SIMD>(toCanonicalInput(x)));
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<bool SIMD>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::CanonicalOutput
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::evalCanonical(const CanonicalInput &x) const noexcept {
    const int nCoeffsRt = (NCOEFFS ? static_cast<int>(NCOEFFS) : static_cast<int>(coeffsMd.extent(0)));
    return poly_eval::horner<NCOEFFS, SIMD, CanonicalOutput>(mapFromDomain(x), coeffsMd, nCoeffsRt);
}
PF_FAST_EVAL_END

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(const CanonicalInput *pts, CanonicalOutput *out,
                                                     std::size_t count) const noexcept {
    for (std::size_t i = 0; i < count; ++i) out[i] = evalCanonical<>(pts[i]);
}

#if defined(__cpp_lib_span) && (__cpp_lib_span >= 202002L)
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(std::span<const CanonicalInput> pts,
                                                     std::span<CanonicalOutput> out) const {
    if (pts.size() != out.size()) {
        throw std::invalid_argument("Input and output spans must have equal length");
    }
    operator()(pts.data(), out.data(), pts.size());
}
#endif

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<class Points, class Outputs, class>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::operator()(const Points &pts, Outputs &out) const {
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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE>::coeffImpl(
    const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffsMd[std::array<std::size_t, DIM + 1>{static_cast<std::size_t>(idx[I])..., k}];
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE>::coeff(
    const IdxArray &idx, std::size_t k) noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<DIM>{});
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<typename IdxArray, std::size_t... I>
constexpr const typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE>::coeffImpl(
    const IdxArray &idx, std::size_t k, std::index_sequence<I...>) const noexcept {
    return coeffsMd[std::array<std::size_t, DIM + 1>{static_cast<std::size_t>(idx[I])..., k}];
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<class IdxArray>
[[nodiscard]] constexpr const typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::Scalar &FuncEvalND<Func, NCOEFFS, FUSION_MODE>::coeff(
    const IdxArray &idx, std::size_t k) const noexcept {
    return coeffImpl<IdxArray>(idx, k, std::make_index_sequence<DIM>{});
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE>::makeExtents(int nCoeffsPerAxis) noexcept -> Extents {
    if constexpr (IS_STATIC) {
        return detail::makeStaticExtents<NCOEFFS, DIM, OUT_DIM>(std::make_index_sequence<DIM>{});
    } else {
        return makeExtents(nCoeffsPerAxis, std::make_index_sequence<DIM + 1>{});
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<std::size_t... Is>
auto FuncEvalND<Func, NCOEFFS, FUSION_MODE>::makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept -> Extents {
    return Extents{(Is < DIM ? static_cast<std::size_t>(nCoeffsPerAxis) : static_cast<std::size_t>(OUT_DIM))...};
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr std::size_t FuncEvalND<Func, NCOEFFS, FUSION_MODE>::storageRequired(const int nCoeffsPerAxis) noexcept {
    auto ext = makeExtents(nCoeffsPerAxis);
    auto mapping = typename Mdspan::mapping_type{ext};
    return mapping.required_span_size();
}
template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::buildCoeffs(int nCoeffsPerAxis, Func f,
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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::convertAxesToMonomialBasis(int nCoeffsPerAxis,
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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::reverseCoefficientOrder(int nCoeffsPerAxis) {
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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<class Point>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::CanonicalInput FuncEvalND<Func, NCOEFFS, FUSION_MODE>::toCanonicalInput(
    const Point &x) noexcept {
    return detail::fixedContainerCast<CanonicalInput>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::InputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::fromCanonicalInput(const CanonicalInput &x) noexcept {
    return detail::fixedContainerCast<InputType>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::CanonicalOutput
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::toCanonicalOutput(const OutputType &x) noexcept {
    return detail::fixedContainerCast<CanonicalOutput>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::OutputType
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::fromCanonicalOutput(const CanonicalOutput &x) noexcept {
    return detail::fixedContainerCast<OutputType>(x);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::CanonicalInput
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::mapToDomain(const CanonicalInput &x, const DomainParams &dp) noexcept {
    if (dp.identityDomain) return x;
    return polyfit::internal::helpers::mapToDomainArray<Scalar, DIM>(x, dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
[[nodiscard]] constexpr typename FuncEvalND<Func, NCOEFFS, FUSION_MODE>::CanonicalInput
FuncEvalND<Func, NCOEFFS, FUSION_MODE>::mapFromDomain(const CanonicalInput &x) const noexcept {
    if constexpr (FUSION_MODE == FusionMode::Always) {
        return x;
    } else {
        if (domain_.identityDomain) return x;
        return polyfit::internal::helpers::mapFromDomainArray<Scalar, DIM>(x, domain_.invSpan, domain_.sumEndpoints);
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
[[nodiscard]] constexpr bool FuncEvalND<Func, NCOEFFS, FUSION_MODE>::shouldFuseAxis(
    const DomainParams &dp, std::size_t axis, int nCoeffsPerAxis) const noexcept {
    if constexpr (FUSION_MODE == FusionMode::Always) {
        return true;
    } else {
        const auto alpha = Scalar(2) * static_cast<Scalar>(dp.invSpan[axis]);
        const auto beta = -static_cast<Scalar>(dp.sumEndpoints[axis]) * static_cast<Scalar>(dp.invSpan[axis]);
        const auto condBase = detail::math::abs(alpha) + detail::math::abs(beta) + Scalar(1);
        constexpr auto maxLog = Scalar(std::numeric_limits<Scalar>::digits10 - 3);
        return nCoeffsPerAxis > 1 &&
               Scalar(nCoeffsPerAxis - 1) * detail::math::log10(condBase) < maxLog;
    }
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::fuseNDDomain(DomainParams &dp, int nCoeffsPerAxis) {
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
        if (!shouldFuseAxis(dp, axis, nCoeffsPerAxis)) continue;

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

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::computeScaling(const InputType &a, const InputType &b,
                                                                      DomainParams &dp) const noexcept {
    polyfit::internal::helpers::computeScalingArray<Scalar, DIM>(toCanonicalInput(a), toCanonicalInput(b), dp.invSpan,
                                                                 dp.sumEndpoints);
    dp.identityDomain = polyfit::internal::helpers::isIdMap(dp.invSpan, dp.sumEndpoints);
}

template<class Func, std::size_t NCOEFFS, FusionMode FUSION_MODE>
template<std::size_t Rank, class F>
constexpr void FuncEvalND<Func, NCOEFFS, FUSION_MODE>::forEachIndex(const std::array<int, Rank> &ext, F &&body) {
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
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported or duplicate fit tag");
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    using Evaluator = FitEvaluator<Func, NCOEFFS, Options::ITERS, Options::FUSION_MODE>;
    return Evaluator(F, a, b);
}

template<class Func, class Spec, class... Tags>
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, Spec spec, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported or duplicate fit tag");

    if constexpr (isIntegralLike_v<Spec>) {
        const auto nCoeffs = detail::validatePositiveCoeffCount(static_cast<int>(spec));
        if constexpr (takesNdInput_v<Func>) {
            return FuncEvalND<Func, 0, Options::FUSION_MODE>(F, nCoeffs, a, b);
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
    static_assert(MAX_NCOEFFS > 0, "Max coefficient count must be positive.");
    static_assert(EVAL_POINTS > 1, "Number of evaluation points must be greater than 1.");

    constexpr auto nCoeffs = [F] {
        constexpr auto computeError = [F](const auto &evaluator) {
            constexpr auto ep = detail::linspace<static_cast<int>(EVAL_POINTS)>(a, b);
            double maxErr = 0.0;
            for (const auto &pt : ep) {
                const auto actual = F(pt);
                const auto approx = evaluator.template operator()<false>(pt);
                maxErr = std::max(detail::relativeL2Norm(approx, actual), maxErr);
            }
            return maxErr;
        };
        int result = 0;
        poet::static_for<1, MAX_NCOEFFS + 1>([&](auto i) {
            if (result != 0) return;
            using Evaluator = std::conditional_t<takesNdInput_v<Func>, FuncEvalND<Func, i, FusionMode::Auto>,
                                                 FuncEval<Func, i, ITERS>>;
            if constexpr (computeError(Evaluator(F, a, b)) <= EPS) {
                result = i;
            }
        });
        return result;
    }();
    static_assert(nCoeffs != 0, "No coefficient count found for requested error tolerance.");
    using Evaluator = std::conditional_t<takesNdInput_v<Func>, FuncEvalND<Func, nCoeffs, FusionMode::Auto>,
                                         FuncEval<Func, nCoeffs, ITERS>>;
    return Evaluator(F, a, b);
}
#endif

template<typename... EvalTypes>
[[nodiscard]] PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...> pack(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval

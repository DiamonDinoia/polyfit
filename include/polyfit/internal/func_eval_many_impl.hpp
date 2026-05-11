#pragma once

#include "impl_common.hpp"

// Out-of-class definitions for poly_eval::FuncEvalMany.

namespace poly_eval {

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
    // Column-tiled SIMD-across-points path only pays off once numPoints
    // amortises blockSize = simdSize * unrollFactor setup and the COUNT-wide
    // re-load of column coefficients. Below the crossover, looping the
    // SIMD-across-polynomials single-point path (horner_transposed via
    // evalInputs) wins.
    constexpr std::size_t simdSize = xsimd::batch<OutputType>::size;
    constexpr std::size_t unrollFactor = detail::optimalManyEvalUf<OutputType>();
    constexpr std::size_t kSmallPointThreshold = simdSize * unrollFactor;

    if (numPoints < kSmallPointThreshold) PF_UNLIKELY {
        for (std::size_t i = 0; i < numPoints; ++i) {
            const auto results = (*this)(x[i]);
            for (std::size_t j = 0; j < COUNT; ++j) out[i * COUNT + j] = results[j];
        }
        return;
    }

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
} // namespace poly_eval

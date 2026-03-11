// Internal domain mapping and coefficient transformation helpers.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "api_types.hpp"
#include "simd_utils.h"

namespace polyfit::internal::helpers {

// Scalar mapping to canonical domain [-1,1]
template<class ArgT, class ScalarT>
constexpr ArgT mapToDomainScalar(const ArgT arg, const ScalarT invSpan, const ScalarT sumEndpoints) noexcept {
    if constexpr (std::is_arithmetic_v<ArgT>) {
        return static_cast<ArgT>(ArgT(0.5) * ((arg / static_cast<ArgT>(invSpan)) + static_cast<ArgT>(sumEndpoints)));
    } else {
        return static_cast<ArgT>(ArgT(0.5) * ((arg / ArgT(invSpan)) + ArgT(sumEndpoints)));
    }
}

// Scalar/batch mapping from canonical domain back to original
template<class ArgT, class ScalarT>
PF_ALWAYS_INLINE constexpr ArgT mapFromDomainScalar(const ArgT arg, const ScalarT invSpan,
                                                    const ScalarT sumEndpoints) noexcept {
    return poly_eval::detail::fma(ArgT(2), arg, -ArgT(sumEndpoints)) * ArgT(invSpan);
}

// Array (std::array) mapping overloads
template<class T, std::size_t N>
constexpr std::array<T, N> mapToDomainArray(const std::array<T, N> &t, const std::array<T, N> &invSpan,
                                            const std::array<T, N> &sumEndpoints) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = mapToDomainScalar<T, T>(t[d], invSpan[d], sumEndpoints[d]);
    }
    return out;
}

template<class T, std::size_t N>
PF_ALWAYS_INLINE constexpr std::array<T, N>
mapFromDomainArray(const std::array<T, N> &x, const std::array<T, N> &invSpan,
                   const std::array<T, N> &sumEndpoints) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = mapFromDomainScalar<T, T>(x[d], invSpan[d], sumEndpoints[d]);
    }
    return out;
}

template<class T, std::size_t N>
constexpr void computeScalingArray(const std::array<T, N> &from, const std::array<T, N> &to,
                                   std::array<T, N> &invSpan, std::array<T, N> &sumEndpoints) noexcept {
    for (std::size_t dim = 0; dim < N; ++dim) {
        invSpan[dim] = T(1) / (to[dim] - from[dim]);
        sumEndpoints[dim] = (to[dim] + from[dim]);
    }
}

namespace eft {

template<class T> constexpr std::pair<T, T> two_sum(T a, T b) noexcept {
    T s = a + b;
    T v = s - a;
    return {s, (a - (s - v)) + (b - v)};
}

template<class T> constexpr std::pair<T, T> two_prod(T a, T b) noexcept {
    T p = a * b;
    PF_IF_CONSTEVAL { return {p, T(0)}; }
    return {p, std::fma(a, b, -p)};
}

} // namespace eft

namespace detail {

template<class T> constexpr void applyTaylorShift(T *coeffs, std::vector<T> &comp, std::size_t n, T beta) noexcept {
    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = n - 2; j >= i; --j) {
            auto [p, ep] = eft::two_prod(beta, coeffs[j + 1]);
            auto [s, es] = eft::two_sum(coeffs[j], p);
            comp[j] += ep + es + beta * comp[j + 1];
            coeffs[j] = s;
            if (j == i) break;
        }
    }
}

template<class T> constexpr void addCompensation(T *coeffs, const std::vector<T> &comp, std::size_t n) noexcept {
    for (std::size_t k = 0; k < n; ++k) coeffs[k] += comp[k];
}

template<class T> constexpr void advancePower(T &alphaPow, T &alphaErr, T alpha) noexcept {
    auto [nextPow, nextErr] = eft::two_prod(alphaPow, alpha);
    alphaErr = alphaErr * alpha + nextErr;
    alphaPow = nextPow;
}

template<class T> constexpr void scaleByPower(T *coeffs, std::size_t n, T alpha) noexcept {
    T alphaPow = alpha;
    T alphaErr = T(0);
    for (std::size_t k = 1; k < n; ++k) {
        auto [p, ep] = eft::two_prod(coeffs[k], alphaPow);
        coeffs[k] = p + (ep + coeffs[k] * alphaErr);
        advancePower(alphaPow, alphaErr, alpha);
    }
}

template<class T>
constexpr T shiftPart(T current, T next, T nextComp, T beta, T &comp) noexcept {
    auto [prod, prodErr] = eft::two_prod(beta, next);
    auto [sum, sumErr] = eft::two_sum(current, prod);
    comp += prodErr + sumErr + beta * nextComp;
    return sum;
}

template<class T>
constexpr T scalePart(T value, T alphaPow, T alphaErr) noexcept {
    auto [prod, prodErr] = eft::two_prod(value, alphaPow);
    return prod + (prodErr + value * alphaErr);
}

template<class T>
constexpr void addCompensation(std::complex<T> *coeffs, const std::vector<T> &compRe, const std::vector<T> &compIm,
                               std::size_t n) noexcept {
    for (std::size_t k = 0; k < n; ++k) coeffs[k] += std::complex<T>(compRe[k], compIm[k]);
}

} // namespace detail

// Fold q(x) = p(alpha * x + beta) into the coefficient array in place.
template<class T, poly_eval::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
constexpr void fuseLinearMap(T *coeffs, std::size_t n, T alpha, T beta) noexcept {
    if (n <= 1) return;

    std::reverse(coeffs, coeffs + n);

    std::vector<T> comp(n, T(0));
    detail::applyTaylorShift(coeffs, comp, n, beta);
    detail::addCompensation(coeffs, comp, n);
    detail::scaleByPower(coeffs, n, alpha);

    std::reverse(coeffs, coeffs + n);
}

template<class T> constexpr void fuseLinearMap(std::complex<T> *coeffs, std::size_t n, T alpha, T beta) noexcept {
    if (n <= 1) return;

    std::reverse(coeffs, coeffs + n);

    std::vector<T> compRe(n, T(0));
    std::vector<T> compIm(n, T(0));

    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = n - 2; j >= i; --j) {
            const auto real = detail::shiftPart(coeffs[j].real(), coeffs[j + 1].real(), compRe[j + 1], beta, compRe[j]);
            const auto imag = detail::shiftPart(coeffs[j].imag(), coeffs[j + 1].imag(), compIm[j + 1], beta, compIm[j]);
            coeffs[j] = std::complex<T>(real, imag);
            if (j == i) break;
        }
    }

    detail::addCompensation(coeffs, compRe, compIm, n);

    T alphaPow = alpha;
    T alphaErr = T(0);
    for (std::size_t k = 1; k < n; ++k) {
        coeffs[k] = std::complex<T>(detail::scalePart(coeffs[k].real(), alphaPow, alphaErr),
                                    detail::scalePart(coeffs[k].imag(), alphaPow, alphaErr));
        detail::advancePower(alphaPow, alphaErr, alpha);
    }

    std::reverse(coeffs, coeffs + n);
}

} // namespace polyfit::internal::helpers

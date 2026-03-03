// Helper utilities for internal use (domain mapping, scaling)
#pragma once

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <cmath>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "simd_utils.h"

namespace polyfit::internal::helpers {

// Scalar mapping to canonical domain [-1,1]
template <class ArgT, class ScalarT>
PF_ALWAYS_INLINE constexpr ArgT map_to_domain_scalar(const ArgT arg, const ScalarT low, const ScalarT hi) noexcept {
    if constexpr (std::is_arithmetic_v<ArgT>) {
        return static_cast<ArgT>(ArgT(0.5) * ((arg / static_cast<ArgT>(low)) + static_cast<ArgT>(hi)));
    }
    return static_cast<ArgT>(ArgT(0.5) * ((arg / ArgT(low)) + ArgT(hi)));
}

// Scalar/batch mapping from canonical domain back to original
template <class ArgT, class ScalarT>
PF_ALWAYS_INLINE constexpr ArgT map_from_domain_scalar(const ArgT arg, const ScalarT low, const ScalarT hi) noexcept {
    return poly_eval::detail::fma(ArgT(2), arg, -ArgT(hi)) * ArgT(low);
}

// Array (std::array) mapping overloads
template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr std::array<T, N>
map_to_domain_array(const std::array<T, N> &t, const std::array<T, N> &low, const std::array<T, N> &hi) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = map_to_domain_scalar<T, T>(t[d], low[d], hi[d]);
    }
    return out;
}

template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr std::array<T, N>
map_from_domain_array(const std::array<T, N> &x, const std::array<T, N> &low, const std::array<T, N> &hi) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = map_from_domain_scalar<T, T>(x[d], low[d], hi[d]);
    }
    return out;
}

template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr void
compute_scaling_array(const std::array<T, N> &from, const std::array<T, N> &to,
                      std::array<T, N> &low, std::array<T, N> &high) noexcept {
    for (std::size_t dim = 0; dim < N; ++dim) {
        low[dim] = T(1) / (to[dim] - from[dim]);
        high[dim] = (to[dim] + from[dim]);
    }
}

//------------------------------------------------------------------------------
// Error-free transformations for compensated arithmetic
//------------------------------------------------------------------------------
namespace eft {

// twoSum: s + e = a + b exactly (Knuth/Moller)
template <class T>
PF_ALWAYS_INLINE constexpr std::pair<T, T> two_sum(T a, T b) noexcept {
    T s = a + b;
    T v = s - a;
    return {s, (a - (s - v)) + (b - v)};
}

// twoProd via FMA: p + e = a * b exactly
template <class T>
PF_ALWAYS_INLINE constexpr std::pair<T, T> two_prod(T a, T b) noexcept {
    T p = a * b;
    return {p, std::fma(a, b, -p)};
}

} // namespace eft

//------------------------------------------------------------------------------
// Fuse linear domain mapping into polynomial coefficients (in-place).
// Given polynomial p(t) in Horner order (highest degree first), computes
// q(x) = p(alpha*x + beta) so that evaluation no longer needs per-point mapping.
//
// Uses compensated (error-free) arithmetic in the Taylor shift and alpha scaling
// to recover machine precision: O(n^2 * eps^2) error instead of O(n^2 * eps).
//------------------------------------------------------------------------------

// Real arithmetic types: full twoProd + twoSum compensation
template <class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
constexpr void fuse_linear_map(T *coeffs, std::size_t n, T alpha, T beta) noexcept {
    if (n <= 1) return;

    std::reverse(coeffs, coeffs + n);

    // --- Compensated Taylor shift: p(x) -> p(x + beta) ---
    // Each step: coeffs[j] += beta * coeffs[j+1]
    // We track per-coefficient compensation using error-free transformations.
    std::vector<T> comp(n, T(0));

    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = n - 2; j >= i; --j) {
            // twoProd: p + ep = beta * coeffs[j+1] exactly
            auto [p, ep] = eft::two_prod(beta, coeffs[j + 1]);
            // First-order correction from compensation of coeffs[j+1]
            T p_comp = beta * comp[j + 1];
            // twoSum: s + es = coeffs[j] + p exactly
            auto [s, es] = eft::two_sum(coeffs[j], p);
            comp[j] += ep + es + p_comp;
            coeffs[j] = s;
            if (j == i) break;
        }
    }

    // Apply accumulated compensation
    for (std::size_t k = 0; k < n; ++k)
        coeffs[k] += comp[k];

    // --- Compensated alpha scaling: coeff[k] *= alpha^k ---
    T alpha_pow = alpha;
    T alpha_err = T(0);
    for (std::size_t k = 1; k < n; ++k) {
        auto [p, ep] = eft::two_prod(coeffs[k], alpha_pow);
        coeffs[k] = p + (ep + coeffs[k] * alpha_err);
        // Accumulate alpha power with error tracking
        auto [ap, ae] = eft::two_prod(alpha_pow, alpha);
        alpha_err = alpha_err * alpha + ae;
        alpha_pow = ap;
    }

    std::reverse(coeffs, coeffs + n);
}

// Complex coefficients with real mapping parameters: component-wise compensation
template <class T>
constexpr void fuse_linear_map(std::complex<T> *coeffs, std::size_t n, T alpha, T beta) noexcept {
    if (n <= 1) return;

    std::reverse(coeffs, coeffs + n);

    // --- Compensated Taylor shift (component-wise) ---
    std::vector<T> comp_re(n, T(0));
    std::vector<T> comp_im(n, T(0));

    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = n - 2; j >= i; --j) {
            // Real part: coeffs[j].real() += beta * coeffs[j+1].real()
            auto [pr, epr] = eft::two_prod(beta, coeffs[j + 1].real());
            T pr_comp = beta * comp_re[j + 1];
            auto [sr, esr] = eft::two_sum(coeffs[j].real(), pr);
            comp_re[j] += epr + esr + pr_comp;

            // Imaginary part: coeffs[j].imag() += beta * coeffs[j+1].imag()
            auto [pi, epi] = eft::two_prod(beta, coeffs[j + 1].imag());
            T pi_comp = beta * comp_im[j + 1];
            auto [si, esi] = eft::two_sum(coeffs[j].imag(), pi);
            comp_im[j] += epi + esi + pi_comp;

            coeffs[j] = std::complex<T>(sr, si);
            if (j == i) break;
        }
    }

    // Apply compensation
    for (std::size_t k = 0; k < n; ++k)
        coeffs[k] += std::complex<T>(comp_re[k], comp_im[k]);

    // --- Compensated alpha scaling (component-wise) ---
    T alpha_pow = alpha;
    T alpha_err = T(0);
    for (std::size_t k = 1; k < n; ++k) {
        auto [pr, epr] = eft::two_prod(coeffs[k].real(), alpha_pow);
        auto [pi, epi] = eft::two_prod(coeffs[k].imag(), alpha_pow);
        coeffs[k] = std::complex<T>(
            pr + (epr + coeffs[k].real() * alpha_err),
            pi + (epi + coeffs[k].imag() * alpha_err));
        auto [ap, ae] = eft::two_prod(alpha_pow, alpha);
        alpha_err = alpha_err * alpha + ae;
        alpha_pow = ap;
    }

    std::reverse(coeffs, coeffs + n);
}

} // namespace polyfit::internal::helpers

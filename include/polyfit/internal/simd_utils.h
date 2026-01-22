#pragma once

// Lightweight shim to centralize xsimd include for internal implementation files.
// Keep this header minimal to avoid duplicating helpers defined in utils.h.

#include <xsimd/xsimd.hpp>
#include <vector>
#include <array>
#include <complex>
#include <type_traits>
#include <cmath>
#include <cstdint>
#include <algorithm>

#include "macros.h"

namespace poly_eval {
namespace detail {

// detect xsimd batches
template <typename T>
struct is_xsimd_batch : std::false_type {};

template <typename T, typename A>
struct is_xsimd_batch<xsimd::batch<T, A>> : std::true_type {};

template <typename Batch>
inline constexpr bool is_xsimd_batch_v = is_xsimd_batch<std::remove_cv_t<Batch>>::value;

 // explicit overload for xsimd::batch to avoid template-deduction ambiguities
template <typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<T, A> fma(const xsimd::batch<T, A> &x,
                                                  const xsimd::batch<T, A> &y,
                                                  const xsimd::batch<T, A> &z) noexcept {
    return xsimd::fma(x, y, z);
}

// mixed batch: complex batch (a) * real batch (b) + complex batch (c)
// supports cases where template deduction would otherwise fail
template <typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(
    const xsimd::batch<std::complex<T>, A> &a,
    const xsimd::batch<T, A> &b,
    const xsimd::batch<std::complex<T>, A> &c) noexcept {
    auto ar = xsimd::real(a);
    auto ai = xsimd::imag(a);
    auto cr = xsimd::real(c);
    auto ci = xsimd::imag(c);
    auto rr = xsimd::fma(ar, b, cr);
    auto ri = xsimd::fma(ai, b, ci);
    return xsimd::batch<std::complex<T>, A>(rr, ri);
}

// mixed batch: complex batch (a) * complex-real batch (b) + complex batch (c)
// (case where b is complex but with real multiplier components)
template <typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(
    const xsimd::batch<std::complex<T>, A> &a,
    const xsimd::batch<std::complex<T>, A> &b,
    const xsimd::batch<std::complex<T>, A> &c) noexcept {
    // delegate to xsimd::fma for complex batches when available
    return xsimd::fma(a, b, c);
}

// arithmetic non-batch types (scalars)
template <typename T, typename = std::enable_if_t<!is_xsimd_batch_v<T> && std::is_arithmetic_v<T>>>
constexpr PF_ALWAYS_INLINE T fma(const T &a, const T &b, const T &c) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
        return std::fma(a, b, c);
    } else {
        return (a * b) + c;
    }
}

// complex<T> (scalar, non-batch)
template <typename T>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const T &b,
                                               const std::complex<T> &c) noexcept {
    return std::complex<T>(std::fma(a.real(), b, c.real()), std::fma(a.imag(), b, c.imag()));
}

// mixed complex / scalar b (scalar, non-batch)
template <typename T, typename B, typename = std::enable_if_t<std::is_arithmetic_v<B>>>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const B &b,
                                               const std::complex<T> &c) noexcept {
    return fma(a, static_cast<T>(b), c);
}

} // namespace detail

// AlignedBuffer: conditional alias that uses xsimd aligned_allocator for dynamic storage.
template <typename T, std::size_t N_compile_time_val, std::size_t alignment>
using AlignedBuffer = std::conditional_t<
    (N_compile_time_val == 0),
    std::vector<T, xsimd::aligned_allocator<T, alignment>>,
    std::array<T, N_compile_time_val>>;

} // namespace poly_eval

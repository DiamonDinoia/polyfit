#pragma once

// Lightweight shim to centralize xsimd include for internal implementation files.
// Keep this header minimal to avoid duplicating helpers defined in utils.h.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <poet/core/register_info.hpp>
#include <type_traits>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "macros.h"

namespace poly_eval {
namespace detail {

// Optimal Horner unroll factor derived from vector register pressure.
// Each lane: 1 pt_batch + 1 acc_batch = 2 vector registers.
// Reserve: 1 broadcast + 2 scratch.
// AVX2/SSE (16 regs) → 6, AVX-512/NEON/SVE (32 regs) → 14.
template<typename T> PF_C23CONSTEVAL std::size_t optimal_horner_uf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    return (nregs - 3) / 2;
}

// Optimal unroll factor for FuncEvalMany bulk eval.
// Matches optimal_horner_uf: UF pt + UF acc + 1 broadcast + 2 scratch = 2*UF + 3.
// Benchmarked: UF=6 beats UF=7 (icache pressure from larger blocks offsets register gain).
// AVX2/SSE (16 vregs) → 6, AVX-512/NEON/SVE (32 vregs) → 14.
template<typename T> PF_C23CONSTEVAL std::size_t optimal_many_eval_uf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    return (nregs - 3) / 2;
}

// detect xsimd batches
template<typename T> struct is_xsimd_batch : std::false_type {};

template<typename T, typename A> struct is_xsimd_batch<xsimd::batch<T, A>> : std::true_type {};

template<typename Batch> inline constexpr bool is_xsimd_batch_v = is_xsimd_batch<std::remove_cv_t<Batch>>::value;

// explicit overload for xsimd::batch to avoid template-deduction ambiguities
template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<T, A> fma(const xsimd::batch<T, A> &x, const xsimd::batch<T, A> &y,
                                                  const xsimd::batch<T, A> &z) noexcept {
    return xsimd::fma(x, y, z);
}

// mixed batch: complex batch (a) * real batch (b) + complex batch (c)
// supports cases where template deduction would otherwise fail
template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(const xsimd::batch<std::complex<T>, A> &a,
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
template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(const xsimd::batch<std::complex<T>, A> &a,
                                                                const xsimd::batch<std::complex<T>, A> &b,
                                                                const xsimd::batch<std::complex<T>, A> &c) noexcept {
    // delegate to xsimd::fma for complex batches when available
    return xsimd::fma(a, b, c);
}

// arithmetic non-batch types (scalars)
// std::fma is constexpr only as a GCC extension in C++20 (standardised in C++26).
// Use a*b+c in constant-evaluated contexts to keep this function constexpr on Clang.
template<typename T, typename = std::enable_if_t<!is_xsimd_batch_v<T> && std::is_arithmetic_v<T>>>
constexpr PF_ALWAYS_INLINE T fma(const T &a, const T &b, const T &c) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
        PF_IF_CONSTEVAL { return a * b + c; } // unfused; acceptable for compile-time polynomial fitting
        return std::fma(a, b, c);
    } else {
        return (a * b) + c;
    }
}

// complex<T> (scalar, non-batch)
template<typename T>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const T &b,
                                               const std::complex<T> &c) noexcept {
    PF_IF_CONSTEVAL { return std::complex<T>(a.real() * b + c.real(), a.imag() * b + c.imag()); }
    return std::complex<T>(std::fma(a.real(), b, c.real()), std::fma(a.imag(), b, c.imag()));
}

// mixed complex / scalar b (scalar, non-batch)
template<typename T, typename B, typename = std::enable_if_t<std::is_arithmetic_v<B>>>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const B &b,
                                               const std::complex<T> &c) noexcept {
    return fma(a, static_cast<T>(b), c);
}

} // namespace detail

// AlignedBuffer: conditional alias that uses xsimd aligned_allocator for dynamic storage.
template<typename T, std::size_t N_compile_time_val, std::size_t alignment>
using AlignedBuffer =
    std::conditional_t<(N_compile_time_val == 0), std::vector<T, xsimd::aligned_allocator<T, alignment>>,
                       std::array<T, N_compile_time_val>>;

} // namespace poly_eval

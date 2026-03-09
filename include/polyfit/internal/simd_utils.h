#pragma once

// Lightweight shim to centralize xsimd include for internal implementation files.
// Keep this header minimal to avoid duplicating helpers defined in utils.h.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <poet/core/cpu_info.hpp>
#include <type_traits>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "macros.h"

namespace poly_eval {
namespace detail {

// Round down x to the nearest multiple of N. Uses bitmask for power-of-2 N,
// otherwise falls back to division (compile-time constant → multiply-shift).
template<std::size_t N>
PF_ALWAYS_INLINE constexpr std::size_t round_down(std::size_t x) noexcept {
    static_assert(N > 0);
    if constexpr ((N & (N - 1)) == 0)
        return x & ~(N - 1);
    else
        return (x / N) * N;
}

// Optimal Horner multi-accumulator unroll factor, tuned by vector register count and width.
// Core Horner loop: UF pt_batches + UF acc_batches + 1 coeff broadcast = 2*UF + 1 regs.
// Wider vectors make spills more expensive (32B for AVX vs 16B for SSE), so we reserve
// more headroom for wider ISAs. Benchmarked UF sweep across SSE/AVX2 × float/double:
//   SSE  (128-bit, 16 regs): UF=8 optimal — 16B spill is cheap    → nregs / 2
//   AVX2 (256-bit, 16 regs): UF=7 optimal — 32B spill hurts       → (nregs - 1) / 2
//   AVX-512 (512-bit, 32 regs): extrapolated UF=15                 → (nregs - 1) / 2
template<typename T> PF_C23CONSTEVAL std::size_t optimal_horner_uf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    constexpr std::size_t vreg_bytes = sizeof(T) * xsimd::batch<T>::size;
    constexpr std::size_t reserved = vreg_bytes <= 16 ? 0 : 1;
    return (nregs - reserved) / 2;
}

// Optimal unroll factor for FuncEvalMany bulk eval (domain mapping + strided scatter).
// Same core budget as optimal_horner_uf, plus extra overhead from domain mapping constants
// (2.0, hi, low broadcasts during load phase) and strided scatter temporaries.
// Benchmarked: the extra per-block work shifts the sweet spot down by ~1 vs optimal_horner_uf.
//   SSE  (128-bit, 16 regs): UF=7    → (nregs - 2) / 2
//   AVX2 (256-bit, 16 regs): UF=6    → (nregs - 3) / 2
//   AVX-512 (512-bit, 32 regs): UF=14 → (nregs - 3) / 2
template<typename T> PF_C23CONSTEVAL std::size_t optimal_many_eval_uf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    constexpr std::size_t vreg_bytes = sizeof(T) * xsimd::batch<T>::size;
    constexpr std::size_t reserved = vreg_bytes <= 16 ? 2 : 3;
    return (nregs - reserved) / 2;
}

// Optimal unroll factor for horner_many scalar interleaving (dynamic_for).
// Each lane runs an independent scalar Horner chain needing 1 FP accumulator.
// Reserve 2 regs for xin + temporaries; divide remaining by 2 for acc + coeff load.
// Cap at 8 to bound dynamic_for binary-tail code size.
//   x86-64   (16 FP regs): UF = min((16-2)/2, 8) = 7
//   AArch64  (32 FP regs): UF = min((32-2)/2, 8) = 8
constexpr std::size_t optimal_horner_many_uf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    constexpr std::size_t uf = (nregs - 2) / 2;
    return uf < 8 ? uf : 8;
}

// horner_transposed scalar: plain for beats dynamic_for (body too light: 1 FMA).


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

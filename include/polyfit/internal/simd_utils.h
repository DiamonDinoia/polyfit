#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <poet/core/cpu_info.hpp>
#include <type_traits>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "api_types.hpp"
#include "macros.h"

namespace poly_eval {
namespace detail {

template<std::size_t N>
PF_ALWAYS_INLINE constexpr std::size_t roundDown(std::size_t x) noexcept {
    static_assert(N > 0);
    if constexpr ((N & (N - 1)) == 0) {
        return x & ~(N - 1);
    } else {
        return (x / N) * N;
    }
}

template<typename T, std::size_t RESERVED> PF_C23CONSTEVAL std::size_t registerLimitedUnroll() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    constexpr std::size_t vreg_bytes = sizeof(T) * xsimd::batch<T>::size;
    constexpr std::size_t actualReserved = vreg_bytes <= 16 ? RESERVED - (RESERVED > 0 ? 1 : 0) : RESERVED;
    return (nregs - actualReserved) / 2;
}

template<typename T> PF_C23CONSTEVAL std::size_t optimalHornerUf() noexcept { return registerLimitedUnroll<T, 1>(); }

template<typename T> PF_C23CONSTEVAL std::size_t optimalManyEvalUf() noexcept { return registerLimitedUnroll<T, 3>(); }

template<class T, std::size_t WIDTH = 1> constexpr std::uint8_t minSimdWidth() {
    if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, WIDTH>>) {
        return minSimdWidth<T, WIDTH * 2>();
    } else {
        return WIDTH;
    }
}

template<class T, std::size_t UPPER, std::size_t WIDTH> constexpr std::size_t optimalSimdWidthImpl() {
    if constexpr (WIDTH * 2 <= UPPER && !std::is_void_v<xsimd::make_sized_batch_t<T, WIDTH * 2>>) {
        return optimalSimdWidthImpl<T, UPPER, WIDTH * 2>();
    } else {
        return WIDTH;
    }
}

template<class T, std::size_t WIDTH_LIMIT> constexpr std::size_t optimalSimdWidth() {
    constexpr std::uint8_t archMax = xsimd::batch<T>::size;
    constexpr std::uint8_t upper = (WIDTH_LIMIT < archMax) ? WIDTH_LIMIT : archMax;
    constexpr std::uint8_t start = minSimdWidth<T>();
    if constexpr (start > upper) {
        return start;
    } else {
        return optimalSimdWidthImpl<T, upper, start>();
    }
}

constexpr std::size_t optimalHornerManyUf() noexcept {
    constexpr std::size_t nregs = poet::vector_register_count();
    constexpr std::size_t uf = (nregs - 2) / 2;
    return uf < 8 ? uf : 8;
}

template<typename T> struct IsXsimdBatch : std::false_type {};

template<typename T, typename A> struct IsXsimdBatch<xsimd::batch<T, A>> : std::true_type {};

template<typename Batch> inline constexpr bool isXsimdBatch_v = IsXsimdBatch<std::remove_cv_t<Batch>>::value;

template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<T, A> fma(const xsimd::batch<T, A> &x, const xsimd::batch<T, A> &y,
                                                  const xsimd::batch<T, A> &z) noexcept {
    return xsimd::fma(x, y, z);
}

template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(const xsimd::batch<std::complex<T>, A> &a,
                                                                const xsimd::batch<T, A> &b,
                                                                const xsimd::batch<std::complex<T>, A> &c) noexcept {
    return xsimd::batch<std::complex<T>, A>(xsimd::fma(xsimd::real(a), b, xsimd::real(c)),
                                            xsimd::fma(xsimd::imag(a), b, xsimd::imag(c)));
}

template<typename T, typename A>
constexpr PF_ALWAYS_INLINE xsimd::batch<std::complex<T>, A> fma(const xsimd::batch<std::complex<T>, A> &a,
                                                                const xsimd::batch<std::complex<T>, A> &b,
                                                                const xsimd::batch<std::complex<T>, A> &c) noexcept {
    return xsimd::fma(a, b, c);
}

template<typename T, typename = enable_if_t<!isXsimdBatch_v<T> && std::is_arithmetic_v<T>>>
constexpr PF_ALWAYS_INLINE T fma(const T &a, const T &b, const T &c) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
        PF_IF_CONSTEVAL { return a * b + c; }
        return std::fma(a, b, c);
    } else {
        return (a * b) + c;
    }
}

template<typename T>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const T &b,
                                               const std::complex<T> &c) noexcept {
    PF_IF_CONSTEVAL { return std::complex<T>(a.real() * b + c.real(), a.imag() * b + c.imag()); }
    return std::complex<T>(std::fma(a.real(), b, c.real()), std::fma(a.imag(), b, c.imag()));
}

template<typename T, typename B, typename = enable_if_t<std::is_arithmetic_v<B>>>
constexpr PF_ALWAYS_INLINE std::complex<T> fma(const std::complex<T> &a, const B &b,
                                               const std::complex<T> &c) noexcept {
    return fma(a, static_cast<T>(b), c);
}

} // namespace detail

template<typename T, std::size_t NCOMPILE, std::size_t ALIGNMENT>
using AlignedBuffer =
    std::conditional_t<NCOMPILE == 0, std::vector<T, xsimd::aligned_allocator<T, ALIGNMENT>>,
                       std::array<T, NCOMPILE>>;

} // namespace poly_eval

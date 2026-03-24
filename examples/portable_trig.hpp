#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>

#include "polyfit/internal/macros.h"
#include "polyfit/internal/numeric_utils.h"
#include <poet/poet.hpp>

#if PF_HAS_CXX20
#include <bit>
#include <numbers>
#endif

namespace polyfit_examples {
namespace portable_trig {

template<typename T> struct default_approx_digits;

template<> struct default_approx_digits<float> : std::integral_constant<int, 7> {};
template<> struct default_approx_digits<double> : std::integral_constant<int, 16> {};
template<> struct default_approx_digits<long double> : std::integral_constant<int, 16> {};
namespace detail {

template<typename T> PF_CXX20_CONSTEXPR T pi() noexcept {
    return static_cast<T>(poly_eval::detail::constants::pi);
}

template<typename T> PF_CXX20_CONSTEXPR T pi_over_2() noexcept { return pi<T>() / static_cast<T>(2); }

template<typename T> PF_CXX20_CONSTEXPR T inv_pi_over_2() noexcept { return static_cast<T>(1) / pi_over_2<T>(); }

inline constexpr std::array<double, 6> sin_coeffs = {
    0x1.5e0f86a545d1fp-33, -0x1.ae6015d2aa2ap-26, 0x1.71de379c19d39p-19,
    -0x1.a01a019e7d0c4p-13, 0x1.1111111110afep-7, -0x1.5555555555554p-3,
};

inline constexpr std::array<double, 6> cos_coeffs = {
    0x1.1c0948cd3683ep-29,  -0x1.27e0fe56c4828p-22, 0x1.a019fc8d5ba89p-16,
    -0x1.6c16c1692bdf3p-10, 0x1.5555555554198p-5,   -0x1.ffffffffffffap-2,
};

// Generated offline on u=t^2 over [0, (pi/4)^2] with the numerically stable
// kernels in examples/portable_trig_generate.cpp. Direct t-domain fits of
// (sin(t) - t) / t^3 and (cos(t) - 1) / t^2 do not reproduce these tables.

template<std::size_t N, std::size_t M>
PF_CXX20_CONSTEVAL std::array<double, N> tail(const std::array<double, M> &src) noexcept {
    static_assert(N <= M, "tail<N> requires N <= M");
    std::array<double, N> out{};
    for (std::size_t i = 0; i < N; ++i) out[i] = src[M - N + i];
    return out;
}

PF_CXX20_CONSTEVAL std::size_t nterms_for_digits(int tol_digits) noexcept {
    return (tol_digits <= 4)  ? 2U
           : (tol_digits <= 6)  ? 3U
           : (tol_digits <= 8)  ? 4U
           : (tol_digits <= 12) ? 5U
                                : 6U;
}

template<int TolDigits> struct trig_coefficients {
    static constexpr auto sin_poly = tail<nterms_for_digits(TolDigits)>(sin_coeffs);
    static constexpr auto cos_poly = tail<nterms_for_digits(TolDigits)>(cos_coeffs);
};

template<typename T> PF_ALWAYS_INLINE T range_reduce_mul_add(T a, T b, T c) noexcept {
#if defined(__FMA__)
    return std::fma(a, b, c);
#else
    return a * b + c;
#endif
}

template<std::size_t N, typename T> PF_ALWAYS_INLINE T eval_poly_horner(const std::array<double, N> &coeffs, T x) noexcept {
    T acc = static_cast<T>(coeffs[0]);
    poet::static_for<1, N>([&](auto i) {
        acc = acc * x + static_cast<T>(coeffs[i]);
    });
    return acc;
}

template<int TolDigits, typename T> PF_ALWAYS_INLINE std::pair<T, T> eval_reduced_sincos(T t) noexcept {
    const T t2 = t * t;
    const T t3 = t2 * t;
    const T sp = eval_poly_horner(trig_coefficients<TolDigits>::sin_poly, t2);
    const T cp = eval_poly_horner(trig_coefficients<TolDigits>::cos_poly, t2);
    return std::pair<T, T>{sp * t3 + t, cp * t2 + static_cast<T>(1)};
}

template<typename To, typename From> PF_ALWAYS_INLINE To bitwise_cast_runtime(const From &src) noexcept {
    static_assert(sizeof(To) == sizeof(From), "bitwise_cast_runtime requires equal-size types");
    static_assert(std::is_trivially_copyable<To>::value, "destination type must be trivially copyable");
    static_assert(std::is_trivially_copyable<From>::value, "source type must be trivially copyable");
#if PF_HAS_CXX20
    return std::bit_cast<To>(src);
#else
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
#endif
}

template<typename T> PF_ALWAYS_INLINE bool is_finite_value(T value) noexcept {
    if constexpr (std::is_same<T, float>::value) {
        constexpr std::uint32_t exponent_mask = 0x7f800000u;
        return (bitwise_cast_runtime<std::uint32_t>(value) & exponent_mask) != exponent_mask;
    } else if constexpr (std::is_same<T, double>::value) {
        constexpr std::uint64_t exponent_mask = 0x7ff0000000000000ull;
        return (bitwise_cast_runtime<std::uint64_t>(value) & exponent_mask) != exponent_mask;
    } else {
        constexpr T max_value = std::numeric_limits<T>::max();
        return value >= -max_value && value <= max_value;
    }
}

template<typename I> PF_ALWAYS_INLINE unsigned quadrant_from_rounded_multiple(I qi) noexcept {
    static_assert(std::is_integral<I>::value, "quadrant_from_rounded_multiple requires an integral type");
    using unsigned_type = std::make_unsigned_t<I>;
    return static_cast<unsigned>(static_cast<unsigned_type>(qi) & unsigned_type{3});
}

template<typename T> PF_ALWAYS_INLINE std::pair<T, T> remap_quadrant(unsigned quadrant, T s1, T c1) noexcept {
    if constexpr (std::is_same<T, float>::value || std::is_same<T, double>::value) {
        using bits_type = std::conditional_t<std::is_same<T, float>::value, std::uint32_t, std::uint64_t>;
        constexpr unsigned sign_shift = std::numeric_limits<bits_type>::digits - 1U;
        const bits_type s_bits = bitwise_cast_runtime<bits_type>(s1);
        const bits_type c_bits = bitwise_cast_runtime<bits_type>(c1);
        const bits_type swap_mask = bits_type(0) - static_cast<bits_type>(quadrant & 1U);
        const bits_type mixed = (s_bits ^ c_bits) & swap_mask;
        const bits_type sin_sign_mask = static_cast<bits_type>(quadrant & 2U) << (sign_shift - 1U);
        const bits_type cos_sign_mask = sin_sign_mask ^ (static_cast<bits_type>(quadrant & 1U) << sign_shift);
        const bits_type sin_bits = (s_bits ^ mixed) ^ sin_sign_mask;
        const bits_type cos_bits = (c_bits ^ mixed) ^ cos_sign_mask;
        return std::pair<T, T>{bitwise_cast_runtime<T>(sin_bits), bitwise_cast_runtime<T>(cos_bits)};
    } else {
        const bool swap = (quadrant & 1U) != 0U;
        const bool sin_neg = (quadrant & 2U) != 0U;
        const bool cos_neg = sin_neg != swap;
        const T sin_mag = swap ? c1 : s1;
        const T cos_mag = swap ? s1 : c1;
        return std::pair<T, T>{sin_neg ? -sin_mag : sin_mag, cos_neg ? -cos_mag : cos_mag};
    }
}

PF_FAST_EVAL_BEGIN
template<int TolDigits> PF_ALWAYS_INLINE std::pair<float, float> sincos_impl(float angle) noexcept {
    if (!is_finite_value(angle)) PF_UNLIKELY {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        return std::pair<float, float>{nan, nan};
    }

    const float qf = std::nearbyint(angle * inv_pi_over_2<float>());
    const std::int32_t qi = static_cast<std::int32_t>(qf);
    float reduced = angle;
    if constexpr (TolDigits >= 6) {
        constexpr float float_pio2_1 = 0x1.920000p+0f;
        constexpr float float_pio2_2 = 0x1.fb4p-12f;
        constexpr float float_pio2_3 = 0x1.4442d2p-24f;
        reduced = range_reduce_mul_add(-qf, float_pio2_1, reduced);
        reduced = range_reduce_mul_add(-qf, float_pio2_2, reduced);
        reduced = range_reduce_mul_add(-qf, float_pio2_3, reduced);
    } else {
        reduced = range_reduce_mul_add(-qf, pi_over_2<float>(), reduced);
    }
    const std::pair<float, float> reduced_sc = eval_reduced_sincos<TolDigits>(reduced);
    return remap_quadrant(quadrant_from_rounded_multiple(qi), reduced_sc.first, reduced_sc.second);
}

template<int TolDigits, typename T> PF_ALWAYS_INLINE std::pair<T, T> sincos_impl(T angle) noexcept {
    static_assert(std::is_floating_point<T>::value, "portable_trig::sincos only supports floating-point inputs");
    if (!is_finite_value(angle)) PF_UNLIKELY {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        return std::pair<T, T>{nan, nan};
    }

    if constexpr (std::is_same<T, double>::value && TolDigits > 14) {
        using index_type = std::int64_t;
        constexpr double medium_inv_pio2 = 0x1.45f306dc9c883p-1;
        constexpr double medium_pio2_1 = 0x1.921fb54400000p+0;
        constexpr double medium_pio2_2 = 0x1.0b4611a626331p-34;
        constexpr double medium_pio2_3 = 0x1.1701b839a2520p-88;
        const double qf = std::nearbyint(angle * medium_inv_pio2);
        const index_type qi = static_cast<index_type>(qf);
        double reduced = angle;
        reduced = range_reduce_mul_add(-qf, medium_pio2_1, reduced);
        reduced = range_reduce_mul_add(-qf, medium_pio2_2, reduced);
        if constexpr (TolDigits > 15) { reduced = range_reduce_mul_add(-qf, medium_pio2_3, reduced); }
        const std::pair<double, double> reduced_sc = eval_reduced_sincos<TolDigits>(reduced);
        return remap_quadrant(quadrant_from_rounded_multiple(qi), reduced_sc.first, reduced_sc.second);
    } else {
        const T qf = std::nearbyint(angle * inv_pi_over_2<T>());
        const std::int64_t qi = static_cast<std::int64_t>(qf);
        const T reduced = range_reduce_mul_add(-qf, pi_over_2<T>(), angle);
        const std::pair<T, T> reduced_sc = eval_reduced_sincos<TolDigits>(reduced);
        return remap_quadrant(quadrant_from_rounded_multiple(qi), reduced_sc.first, reduced_sc.second);
    }
}
PF_FAST_EVAL_END

} // namespace detail

template<typename T, int TolDigits = default_approx_digits<T>::value>
PF_ALWAYS_INLINE std::pair<T, T> sincos(T angle) noexcept {
    return detail::sincos_impl<TolDigits>(angle);
}

template<typename T, int TolDigits = default_approx_digits<T>::value>
PF_ALWAYS_INLINE std::complex<T> cis(T angle, T magnitude = T(1)) noexcept {
    const std::pair<T, T> sc = sincos<T, TolDigits>(angle);
    return std::complex<T>(magnitude * sc.second, magnitude * sc.first);
}

template<typename T, int TolDigits = default_approx_digits<T>::value>
PF_ALWAYS_INLINE std::complex<T> polar(T magnitude, T angle) noexcept {
    return cis<T, TolDigits>(angle, magnitude);
}

} // namespace portable_trig
} // namespace polyfit_examples

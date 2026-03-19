#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>
#include <utility>

#include "polyfit/internal/macros.h"
#include "polyfit/internal/numeric_utils.h"

#include <poet/poet.hpp>

#if PF_HAS_CXX20
#include <numbers>
#endif

namespace polyfit_examples {
namespace portable_trig {

template<typename T> struct default_approx_digits;

template<> struct default_approx_digits<float> : std::integral_constant<int, 7> {};
template<> struct default_approx_digits<double> : std::integral_constant<int, 15> {};
template<> struct default_approx_digits<long double> : std::integral_constant<int, 15> {};

namespace detail {

#if defined(__FMA__)
inline constexpr bool kPreferHardwareFma = true;
#else
inline constexpr bool kPreferHardwareFma = false;
#endif

template<typename T> PF_CXX20_CONSTEXPR T pi() noexcept {
#if PF_HAS_CXX20
    return std::numbers::pi_v<T>;
#else
    return static_cast<T>(poly_eval::detail::constants::pi);
#endif
}

template<typename T> PF_CXX20_CONSTEXPR T pi_over_2() noexcept { return pi<T>() / static_cast<T>(2); }

template<typename T> PF_CXX20_CONSTEXPR T inv_pi_over_2() noexcept { return static_cast<T>(1) / pi_over_2<T>(); }

inline constexpr std::array<double, 6> sin_coeffs = {
    0x1.5e585f68f956ep-33, -0x1.ae5f687b275b3p-26, 0x1.71de33799ebc6p-19,
    -0x1.a01a019367fdp-13, 0x1.1111111104f1dp-7,   -0x1.555555555541bp-3,
};

inline constexpr std::array<double, 6> cos_coeffs = {
    0x1.1b88ad1c62723p-29,  -0x1.27df3a1e26a95p-22, 0x1.a019f7fcecefp-16,
    -0x1.6c16c163eaf27p-10, 0x1.555555554ef27p-5,   -0x1.fffffffffff91p-2,
};

// These tables were generated offline with polyfit on the reduced interval
// [a, b] = [-pi/4, pi/4], then rewritten into the
//   sin(t) ~= t + t^3 * P(t^2)
//   cos(t) ~= 1 + t^2 * Q(t^2)
// form used below. A representative generation sketch is:
//
//   using poly_eval::FuseAlways;
//   const double a = -poly_eval::detail::constants::pi / 4.0;
//   const double b =  poly_eval::detail::constants::pi / 4.0;
//   const auto sin_tail = poly_eval::fit<6>(
//       [](double t) { return (std::sin(t) - t) / (t * t * t); },
//       a, b, FuseAlways{});
//   const auto cos_tail = poly_eval::fit<6>(
//       [](double t) { return (std::cos(t) - 1.0) / (t * t); },
//       a, b, FuseAlways{});
//
// The stored coefficients are the Horner-order coefficients of P and Q in u=t^2.
//
// Local performance note from this example on one x86-64 machine, compared to
// std::sin(x) + std::cos(x) in the same harness:
//   -march=x86-64    : roughly 1.2x to 1.6x faster
//   -march=x86-64-v2 : roughly 1.6x to 1.75x faster
//   -march=x86-64-v3 : roughly 1.55x to 1.8x faster
// These numbers are workload- and machine-dependent; rerun portable_trig.cpp on
// the target machine if you need a trustworthy comparison.

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

template<int TolDigits> struct approx_coefficients {
    static constexpr auto sin_poly = tail<nterms_for_digits(TolDigits)>(sin_coeffs);
    static constexpr auto cos_poly = tail<nterms_for_digits(TolDigits)>(cos_coeffs);
};

template<typename T> PF_ALWAYS_INLINE T range_reduce_mul_add(T a, T b, T c) noexcept {
    if constexpr (kPreferHardwareFma && (std::is_same<T, float>::value || std::is_same<T, double>::value)) {
        return poly_eval::detail::math::fma(a, b, c);
    } else {
        return a * b + c;
    }
}

template<typename T> PF_ALWAYS_INLINE T plain_mul_add(T a, T b, T c) noexcept {
    return a * b + c;
}

template<std::size_t N, typename T> PF_ALWAYS_INLINE T eval_horner_baseline(const std::array<double, N> &coeffs, T x) noexcept {
    T acc = static_cast<T>(coeffs[0]);
    poet::static_for<1, N>([&](auto i) {
        acc = plain_mul_add(acc, x, static_cast<T>(coeffs[i]));
    });
    return acc;
}

template<int TolDigits, typename T> PF_ALWAYS_INLINE std::pair<T, T> evaluate_reduced_baseline(T t) noexcept {
    const T t2 = t * t;
    const T t3 = t2 * t;
    const T sp = eval_horner_baseline(approx_coefficients<TolDigits>::sin_poly, t2);
    const T cp = eval_horner_baseline(approx_coefficients<TolDigits>::cos_poly, t2);
    return std::pair<T, T>{plain_mul_add(sp, t3, t), plain_mul_add(cp, t2, static_cast<T>(1))};
}

template<typename T> PF_ALWAYS_INLINE std::pair<T, T> nan_pair() noexcept {
    const T nan = std::numeric_limits<T>::quiet_NaN();
    return std::pair<T, T>{nan, nan};
}

PF_FAST_EVAL_BEGIN
template<int TolDigits, typename T> PF_ALWAYS_INLINE std::pair<T, T> sincos_impl(T angle) noexcept {
    static_assert(std::is_floating_point<T>::value, "portable_trig::sincos only supports floating-point inputs");

    if (!std::isfinite(angle)) return nan_pair<T>();

    using work_type = typename std::conditional<std::is_same<T, float>::value, double, T>::type;
    const work_type angle_w = static_cast<work_type>(angle);
    const work_type qf = std::nearbyint(angle_w * inv_pi_over_2<work_type>());
    const long long qi = static_cast<long long>(qf);
    const work_type reduced = range_reduce_mul_add(-qf, pi_over_2<work_type>(), angle_w);
    const std::pair<work_type, work_type> reduced_sc = evaluate_reduced_baseline<TolDigits>(reduced);
    const T s1 = static_cast<T>(reduced_sc.first);
    const T c1 = static_cast<T>(reduced_sc.second);

    int quadrant = static_cast<int>(qi % 4LL);
    if (quadrant < 0) quadrant += 4;

    if (quadrant == 0) return std::pair<T, T>{s1, c1};
    if (quadrant == 1) return std::pair<T, T>{c1, -s1};
    if (quadrant == 2) return std::pair<T, T>{-s1, -c1};
    return std::pair<T, T>{-c1, s1};
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

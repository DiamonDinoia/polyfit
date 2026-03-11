#pragma once

#if defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L)
#include <bit>
#endif

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "macros.h"

namespace poly_eval::detail {

template<typename T> constexpr auto countr_zero(T x) noexcept {
    static_assert(std::is_unsigned_v<T>, "countr_zero requires an unsigned integral type");
#if defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L)
    return std::countr_zero(x);
#else
    constexpr int width = std::numeric_limits<T>::digits;
    if (x == 0) {
        return width;
    }
    int count = 0;
    while ((x & T{1}) == T{0}) {
        x >>= 1;
        ++count;
    }
    return count;
#endif
}

template<typename T> constexpr std::size_t getAlignment(const T *ptr) noexcept {
    const auto address = reinterpret_cast<std::uintptr_t>(ptr);
    if (address == 0) {
        return 0;
    }
    return static_cast<std::size_t>(1) << countr_zero(address);
}

namespace constants {
inline constexpr double pi = 3.14159265358979323846;
} // namespace constants

namespace math {

template<typename T> constexpr T fma(T a, T b, T c) noexcept {
    PF_IF_CONSTEVAL { return a * b + c; }
    return std::fma(a, b, c);
}

template<typename T> constexpr T abs(T x) noexcept {
    PF_IF_CONSTEVAL { return x < T(0) ? -x : x; }
    return std::abs(x);
}

template<typename T> constexpr T log10(T x) noexcept {
    PF_IF_CONSTEVAL {
        if (x <= T(0)) return T(0);
        T result = T(0);
        while (x >= T(10)) {
            x /= T(10);
            result += T(1);
        }
        while (x < T(1)) {
            x *= T(10);
            result -= T(1);
        }
        const T y = (x - T(1)) / (x + T(1));
        const T y2 = y * y;
        const T lnX = T(2) * y *
                      (T(1) + y2 * (T(1.0 / 3.0) +
                                    y2 * (T(1.0 / 5.0) + y2 * (T(1.0 / 7.0) + y2 * (T(1.0 / 9.0) + y2 / T(11))))));
        return result + lnX / T(2.302585092994046);
    }
    return std::log10(x);
}

constexpr double sqrt(double x) noexcept {
    PF_IF_CONSTEVAL {
        if (x <= 0.0) return 0.0;
        double estimate = x;
        for (int i = 0; i < 100; ++i) {
            const double next = 0.5 * (estimate + x / estimate);
            if (next == estimate) break;
            estimate = next;
        }
        return estimate;
    }
    return std::sqrt(x);
}

} // namespace math

constexpr double cos(const double x) noexcept {
    constexpr double PIO2_HI = 1.57079632679489655800e+00;
    constexpr double PIO2_LO = 6.12323399573676603587e-17;
    constexpr double INV_PIO2 = 6.36619772367581382433e-01;

    const double fn = x * INV_PIO2;
    const int n = static_cast<int>(fn + (fn >= 0.0 ? 0.5 : -0.5));
    const int quadrant = n & 3;
    const double yn = static_cast<double>(-n);
    const double y = (yn * PIO2_HI + x) + yn * PIO2_LO;

    const auto cosPoly = [](const double yy) noexcept {
        constexpr double c1 = 4.16666666666666019037e-02;
        constexpr double c2 = -1.38888888888741095749e-03;
        constexpr double c3 = 2.48015872894767294178e-05;
        constexpr double c4 = -2.75573143513906633035e-07;
        constexpr double c5 = 2.08757232129817482790e-09;
        constexpr double c6 = -1.13596475577881948265e-11;
        const double z = yy * yy;
        double r = c6 * z + c5;
        r = r * z + c4;
        r = r * z + c3;
        r = r * z + c2;
        r = r * z + c1;
        return z * z * r + (1.0 - 0.5 * z);
    };

    const auto sinPoly = [](const double yy) noexcept {
        constexpr double s1 = -1.66666666666666307295e-01;
        constexpr double s2 = 8.33333333332211858878e-03;
        constexpr double s3 = -1.98412698295895385996e-04;
        constexpr double s4 = 2.75573136213857245213e-06;
        constexpr double s5 = -2.50507477628578072866e-08;
        constexpr double s6 = 1.58962301576546568060e-10;
        const double z = yy * yy;
        double r = s6 * z + s5;
        r = r * z + s4;
        r = r * z + s3;
        r = r * z + s2;
        r = r * z + s1;
        return yy * z * r + yy;
    };

    switch (quadrant) {
    case 0:
        return cosPoly(y);
    case 1:
        return -sinPoly(y);
    case 2:
        return -cosPoly(y);
    default:
        return sinPoly(y);
    }
}

} // namespace poly_eval::detail

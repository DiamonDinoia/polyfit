#pragma once
#if defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L)
#include <bit>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <xsimd/xsimd.hpp>

#include "helpers.h"
#include "macros.h"

#if __cplusplus < 202002L
namespace std {
template<typename T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
} // namespace std
#endif

namespace poly_eval {

template<class Func, std::size_t, std::size_t> class FuncEval;
// -----------------------------------------------------------------------------
// Buffer: Conditional type alias for std::vector or std::array
// -----------------------------------------------------------------------------
template<typename T, std::size_t N_compile_time_val>
using Buffer = std::conditional_t<N_compile_time_val == 0, std::vector<T>, std::array<T, N_compile_time_val>>;

// Create a Buffer<T,N>, resizing to runtime_size when N==0 (dynamic).
template<typename T, std::size_t N> constexpr Buffer<T, N> make_buffer(std::size_t runtime_size) {
    Buffer<T, N> buf{};
    if constexpr (N == 0) buf.resize(runtime_size);
    return buf;
}

// -----------------------------------------------------------------------------
// function_traits: Helper to deduce input and output types from a callable
// -----------------------------------------------------------------------------
template<typename T> struct function_traits : function_traits<decltype(&std::remove_cvref_t<T>::operator())> {};

template<typename R, typename Arg> struct function_traits<R (*)(Arg)> {
    using result_type = R;
    using arg0_type = Arg;
};

// std::function<R(Arg)>
template<typename R, typename Arg> struct function_traits<std::function<R(Arg)>> {
    using result_type = R;
    using arg0_type = Arg;
};

// pointer to const member (for lambdas with const operator())
template<typename F, typename R, typename Arg> struct function_traits<R (F::*)(Arg) const> {
    using result_type = R;
    using arg0_type = Arg;
};

// pointer to non‐const member (if you ever need it)
template<typename F, typename R, typename Arg> struct function_traits<R (F::*)(Arg)> {
    using result_type = R;
    using arg0_type = Arg;
};

template<typename R, typename T> struct function_traits<R (*)(const T &)> {
    using result_type = R;
    using arg0_type = T;
};

template<typename T, typename = void> struct is_tuple_like : std::false_type {};

template<typename T>
struct is_tuple_like<T, std::void_t<decltype(std::tuple_size_v<std::remove_cvref_t<T>>)>> : std::true_type {};

#if __cpp_concepts >= 201907L
template<typename T>
concept tuple_like = is_tuple_like<T>::value;
#endif

// Convenience: size-or-zero that never hard-errors
template<typename T, typename = void> struct tuple_size_or_zero : std::integral_constant<std::size_t, 0> {};

template<typename T>
struct tuple_size_or_zero<T, std::void_t<decltype(std::tuple_size_v<std::remove_cvref_t<T>>)>>
    : std::integral_constant<std::size_t, std::tuple_size_v<std::remove_cvref_t<T>>> {};

template<typename T> struct is_func_eval : std::false_type {};
template<typename... Args> struct is_func_eval<FuncEval<Args...>> : std::true_type {};
} // namespace poly_eval

template<typename T> struct is_complex : std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : std::true_type {};
template<typename T> inline constexpr bool is_complex_v = is_complex<std::remove_cv_t<T>>::value;

// —— detect whether T has a std::tuple_size<T>::value member (e.g. std::array, tuple, etc.)
template<typename, typename = void> struct has_tuple_size : std::false_type {};
template<typename T> struct has_tuple_size<T, std::void_t<decltype(std::tuple_size<T>::value)>> : std::true_type {};
template<typename T> inline constexpr bool has_tuple_size_v = has_tuple_size<std::remove_cv_t<T>>::value;

// Safely get value_type for containers, or return T for scalars.
template<typename T, typename = void> struct value_type_or_identity {
    using type = T;
};

template<typename T> struct value_type_or_identity<T, std::void_t<typename T::value_type>> {
    using type = typename T::value_type;
};

namespace poly_eval::detail {

namespace eft = polyfit::internal::helpers::eft;

// std::countr_zero returns the number of trailing zero bits.
// If an address is N-byte aligned, its N lowest bits must be zero.
// So, if an address is 8-byte aligned (e.g., 0x...1000), it has 3 trailing zeros.
// 2^3 = 8.
template<typename T> constexpr auto countr_zero(T x) noexcept {
    static_assert(std::is_unsigned_v<T>, "cntz requires an unsigned integral type");
    static_assert(std::is_unsigned<T>::value, "countr_zero_impl requires an unsigned type");
#if defined(__cpp_lib_bitops) && (__cpp_lib_bitops >= 201907L)
    // C++20: hand work to the standard library
    return std::countr_zero(x);
#else
    // constexpr portable fallback
    constexpr int w = std::numeric_limits<T>::digits;
    if (x == 0) {
        return w;
    }
    int n = 0;
    while ((x & T{1}) == T{0}) {
        x >>= 1;
        ++n;
    }
    return n;
#endif
}

template<typename T> constexpr size_t get_alignment(const T *ptr) noexcept {
    const auto address = reinterpret_cast<std::uintptr_t>(ptr);
    if (address == 0) {
        // A null pointer (or an address of 0) doesn't have a meaningful alignment
        // in the context of data access.
        return 0;
    }
    return static_cast<size_t>(1) << detail::countr_zero(address);
}

namespace constants {
inline constexpr double pi = 3.14159265358979323846;
} // namespace constants

// constexpr-safe fma: std::fma is constexpr only as a GCC extension (standardised in C++26).
// In constant-evaluated contexts fall back to a*b+c to keep PF_C20CONSTEXPR functions
// evaluable at compile time on Clang.
template<typename T> constexpr T ce_fma(T a, T b, T c) noexcept {
    if (PF_IS_CONSTANT_EVALUATED()) return a * b + c;
    return std::fma(a, b, c);
}

constexpr double cos(const double x) noexcept {
    // Constexpr Cody-Waite cos approximation.
    // Uses a*b+c instead of std::fma for MSVC constexpr compatibility.
    // NOTE: do NOT dispatch to std::cos at runtime in C++20. The Cody-Waite
    // values produce consistent Chebyshev nodes that avoid near-cancellation in
    // the Björck-Pereyra divided differences at high degrees (≥ 44). std::cos
    // returns values that differ by up to 1 ULP, which is enough to cause
    // accuracy degradation at high polynomial degrees.
    constexpr double PIO2_HI = 1.57079632679489655800e+00;
    constexpr double PIO2_LO = 6.12323399573676603587e-17;
    constexpr double INV_PIO2 = 6.36619772367581382433e-01;

    const double fn = x * INV_PIO2;
    const int n = static_cast<int>(fn + (fn >= 0.0 ? 0.5 : -0.5));
    const int q = n & 3; // quadrant 0‥3
    const double yn = static_cast<double>(-n);
    const double y = (yn * PIO2_HI + x) + yn * PIO2_LO;

    /* cos & sin minimax polynomials as lambdas with embedded coeffs */
    const auto cos_poly = [](const double yy) noexcept {
        /*
         * The coefficients c1-c6 are under the following license:
         * ====================================================
         * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
         *
         * Developed at SunSoft, a Sun Microsystems, Inc. business.
         * Permission to use, copy, modify, and distribute this
         * software is freely granted, provided that this notice
         * is preserved.
         * ====================================================
         */
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

    const auto sin_poly = [](const double yy) noexcept {
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

    /* quadrant dispatch—only compute what we need */
    switch (q) {
    case 0:
        return cos_poly(y);
    case 1:
        return -sin_poly(y);
    case 2:
        return -cos_poly(y);
    default:
        return sin_poly(y);
    }
}

// GCC false positive: inlining vector alloc/dealloc through bjorck_pereyra and
// newton_to_monomial triggers -Wfree-nonheap-object / -Wnull-dereference.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#endif
// stand-alone Bjorck–Pereyra divided‐difference
template<std::size_t N, class X, class Y>
PF_C20CONSTEXPR Buffer<Y, N> bjorck_pereyra(const Buffer<X, N> &x, const Buffer<Y, N> &y) {
    const std::size_t n = (N == 0 ? x.size() : N);
    Buffer<Y, N> a = y;

    if constexpr (std::is_arithmetic_v<Y>) {
        // Compensated path for real types
        auto comp = make_buffer<Y, N>(n);
        for (auto &v : comp) v = Y(0);

        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                auto [s, es] = eft::two_sum(a[i], -a[i - 1]);
                Y sub_comp = es + (comp[i] - comp[i - 1]);
                Y d = Y(x[i] - x[i - k - 1]);
                Y q = s / d;
                Y r = ce_fma(-q, d, s); // division residual
                comp[i] = (r + sub_comp) / d;
                a[i] = q;
            }
            // Fold compensation back after each outer iteration (lossless)
            for (std::size_t i = n - 1; i > k; --i) {
                auto [s, e] = eft::two_sum(a[i], comp[i]);
                a[i] = s;
                comp[i] = e;
            }
        }
    } else if constexpr (is_complex_v<Y>) {
        // Compensated path for complex types
        // Divisor d is real (type X), so division is component-wise
        using Scalar = typename Y::value_type;
        auto comp_re = make_buffer<Scalar, N>(n);
        auto comp_im = make_buffer<Scalar, N>(n);
        for (auto &v : comp_re) v = Scalar(0);
        for (auto &v : comp_im) v = Scalar(0);

        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                Scalar d = Scalar(x[i] - x[i - k - 1]);
                // Real component
                auto [sr, esr] = eft::two_sum(a[i].real(), -a[i - 1].real());
                Scalar sub_comp_re = esr + (comp_re[i] - comp_re[i - 1]);
                Scalar qr = sr / d;
                Scalar rr = ce_fma(-qr, d, sr);
                comp_re[i] = (rr + sub_comp_re) / d;
                // Imaginary component
                auto [si, esi] = eft::two_sum(a[i].imag(), -a[i - 1].imag());
                Scalar sub_comp_im = esi + (comp_im[i] - comp_im[i - 1]);
                Scalar qi = si / d;
                Scalar ri = ce_fma(-qi, d, si);
                comp_im[i] = (ri + sub_comp_im) / d;

                a[i] = Y(qr, qi);
            }
            // Fold compensation back after each outer iteration (lossless)
            for (std::size_t i = n - 1; i > k; --i) {
                auto [sr, er] = eft::two_sum(a[i].real(), comp_re[i]);
                auto [si, ei] = eft::two_sum(a[i].imag(), comp_im[i]);
                a[i] = Y(sr, si);
                comp_re[i] = er;
                comp_im[i] = ei;
            }
        }
    } else {
        // Fallback: uncompensated
        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                a[i] = (a[i] - a[i - 1]) / static_cast<Y>(x[i] - x[i - k - 1]);
            }
        }
    }

    return a;
}

// stand-alone Newton→monomial conversion
// The inner loop touches c[deg] where deg can reach n, so the workspace
// needs n+1 entries; only c[0..n-1] are meaningful in the result.
template<std::size_t N, class X, class Y>
PF_C20CONSTEXPR Buffer<Y, N> newton_to_monomial(const Buffer<Y, N> &alpha, const Buffer<X, N> &nodes) {
    const std::size_t n = alpha.size();
    constexpr std::size_t WN = (N == 0) ? 0 : N + 1;
    auto c = make_buffer<Y, WN>(n + 1);
    for (auto &v : c) v = Y(0);

    if constexpr (std::is_arithmetic_v<Y>) {
        // Compensated path for real types
        auto comp = make_buffer<Y, WN>(n + 1);
        for (auto &v : comp) v = Y(0);

        std::size_t deg = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++deg;
            for (std::size_t j = deg; j >= 1; --j) {
                auto [p, ep] = eft::two_prod(Y(nodes[ii]), c[j]);
                Y p_comp = Y(nodes[ii]) * comp[j];
                auto [s, es] = eft::two_sum(c[j - 1], -p);
                comp[j] = comp[j - 1] + es - ep - p_comp;
                c[j] = s;
            }
            // c[0] = (-nodes[ii]*c[0]) + alpha[ii]
            auto [p0, ep0] = eft::two_prod(Y(nodes[ii]), c[0]);
            Y p0_comp = Y(nodes[ii]) * comp[0];
            auto [s0, es0] = eft::two_sum(alpha[ii], -p0);
            comp[0] = es0 - ep0 - p0_comp;
            c[0] = s0;
        }
        // Apply compensation
        for (std::size_t k = 0; k <= n; ++k) c[k] += comp[k];
    } else if constexpr (is_complex_v<Y>) {
        // Compensated path for complex types (component-wise)
        using Scalar = typename Y::value_type;
        auto comp_re = make_buffer<Scalar, WN>(n + 1);
        auto comp_im = make_buffer<Scalar, WN>(n + 1);
        for (auto &v : comp_re) v = Scalar(0);
        for (auto &v : comp_im) v = Scalar(0);

        std::size_t deg = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++deg;
            Scalar node = Scalar(nodes[ii]);
            for (std::size_t j = deg; j >= 1; --j) {
                // Real component
                auto [pr, epr] = eft::two_prod(node, c[j].real());
                Scalar pr_comp = node * comp_re[j];
                auto [sr, esr] = eft::two_sum(c[j - 1].real(), -pr);
                comp_re[j] = comp_re[j - 1] + esr - epr - pr_comp;
                // Imaginary component
                auto [pi, epi] = eft::two_prod(node, c[j].imag());
                Scalar pi_comp = node * comp_im[j];
                auto [si, esi] = eft::two_sum(c[j - 1].imag(), -pi);
                comp_im[j] = comp_im[j - 1] + esi - epi - pi_comp;
                c[j] = Y(sr, si);
            }
            // c[0] = (-nodes[ii]*c[0]) + alpha[ii]
            auto [pr0, epr0] = eft::two_prod(node, c[0].real());
            Scalar pr0_comp = node * comp_re[0];
            auto [sr0, esr0] = eft::two_sum(alpha[ii].real(), -pr0);
            comp_re[0] = esr0 - epr0 - pr0_comp;

            auto [pi0, epi0] = eft::two_prod(node, c[0].imag());
            Scalar pi0_comp = node * comp_im[0];
            auto [si0, esi0] = eft::two_sum(alpha[ii].imag(), -pi0);
            comp_im[0] = esi0 - epi0 - pi0_comp;

            c[0] = Y(sr0, si0);
        }
        // Apply compensation
        for (std::size_t k = 0; k <= n; ++k) c[k] += Y(comp_re[k], comp_im[k]);
    } else {
        // Fallback: uncompensated
        std::size_t deg = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++deg;
            for (std::size_t j = deg; j >= 1; --j) c[j] = c[j - 1] - (nodes[ii] * c[j]);
            c[0] = (-nodes[ii] * c[0]) + alpha[ii];
        }
    }

    if constexpr (N == 0) {
        c.resize(n);
        return c;
    } else {
        Buffer<Y, N> result{};
        std::copy_n(c.begin(), N, result.begin());
        return result;
    }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

template<class T, std::size_t N = 1> constexpr std::uint8_t min_simd_width() {
    if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return min_simd_width<T, N * 2>();
    } else {
        return N;
    }
}

template<class T, std::size_t Upper, std::size_t Width> constexpr std::size_t optimal_impl() {
    if constexpr (Width * 2 <= Upper && !std::is_void_v<xsimd::make_sized_batch_t<T, Width * 2>>) {
        return optimal_impl<T, Upper, Width * 2>();
    } else {
        return Width;
    }
}

template<class T, std::size_t N> constexpr std::size_t optimal_simd_width() {
    constexpr std::uint8_t arch_max = xsimd::batch<T>::size;
    constexpr std::uint8_t upper = (N < arch_max) ? N : arch_max;
    constexpr std::uint8_t start = min_simd_width<T>();
    if constexpr (start > upper) {
        // N is smaller than the smallest SIMD batch → scalar
        return start;
    } else {
        return optimal_impl<T, upper, start>();
    }
}

// --- Helper to create static extents when N_compile > 0 --------------------
// 1) free helper: build a static extents type for N_compile > 0
template<std::size_t N_compile, std::size_t Dim, std::size_t OutDim, std::size_t... Is>
PF_C23CONSTEVAL auto make_static_extents_impl(std::index_sequence<Is...>) {
    // expands to extents<N_compile, N_compile, …, OutDim>
    return stdex::extents<std::size_t, ((void)Is, N_compile)..., OutDim>{};
}

template<std::size_t N_compile, std::size_t Dim, std::size_t OutDim>
using static_extents_t = decltype(make_static_extents_impl<N_compile, Dim, OutDim>(std::make_index_sequence<Dim>{}));

// 2) free helper: compute how many entries that layout_right mdspan needs
template<class Scalar, std::size_t N_compile, std::size_t Dim, std::size_t OutDim>
PF_C23CONSTEVAL std::size_t storage_required() {
    // pick the extents type directly (no const!)
    using extents_t = static_extents_t<N_compile, Dim, OutDim>;
    using mdspan_t = stdex::mdspan<Scalar, extents_t, stdex::layout_right>;
    // default‑construct an extents_t and query its mapping
    return typename mdspan_t::mapping_type{extents_t{}}.required_span_size();
}

template<std::size_t N_compile, std::size_t DimIn, std::size_t DimOut, std::size_t... Is>
PF_C23CONSTEVAL auto make_static_extents(std::index_sequence<Is...>) {
    return stdex::extents<std::size_t, ((void)Is, N_compile)..., DimOut>{};
}

template<std::size_t M = 0, typename T>
PF_C20CONSTEXPR auto linspace(const T &start, const T &end, int num_points = M) {
    // we'll store each “point” in a Buffer<T,M>
    Buffer<T, M> points{};
    if (num_points <= 0) {
        return points;
    }

    const auto count = static_cast<std::size_t>(num_points);
    if constexpr (M == 0) {
        points.resize(count);
    }

    // scalar case
    if constexpr (std::is_arithmetic_v<T>) {
        if (num_points <= 1) {
            if (num_points == 1) {
                points[0] = start;
            }
            return points;
        }
        const T step = (end - start) / static_cast<T>(count - 1);
        for (std::size_t i = 0; i < count; ++i) {
            points[i] = start + (static_cast<T>(i) * step);
        }
        return points;
    }
    // array case: T must be std::array<Scalar,D>
    else {
        constexpr std::size_t D = std::tuple_size_v<std::remove_cvref_t<T>>;
        if (num_points <= 1) {
            if (num_points == 1) {
                points[0] = start;
            }
            return points;
        }
        using Scalar = std::remove_cv_t<decltype(start[0])>;
        T step{};
        for (std::size_t i = 0; i < D; ++i) {
            step[i] = (end[i] - start[i]) / static_cast<Scalar>(count - 1);
        }

        for (std::size_t k = 0; k < count; ++k) {
            for (std::size_t i = 0; i < D; ++i) {
                points[k][i] = start[i] + (static_cast<Scalar>(k) * step[i]);
            }
        }
        return points;
    }
}

// constexpr sqrt via Newton-Raphson: std::sqrt is not constexpr until C++26.
constexpr double ce_sqrt(double x) noexcept {
    if (x <= 0.0) return 0.0;
    double r = x;
    for (int i = 0; i < 100; ++i) {
        const double next = 0.5 * (r + x / r);
        if (next == r) break;
        r = next;
    }
    return r;
}

template<typename T> PF_C20CONSTEXPR double relative_error(const T &approx, const T &actual) {
    // std::abs for floats is not constexpr until C++23; use manual abs in constexpr path.
    const auto ce_abs = [](double v) constexpr noexcept {
        return v < 0.0 ? -v : v;
    };
    if constexpr (has_tuple_size_v<T>) {
        double err = 0.0;
        for (std::size_t i = 0; i < std::tuple_size_v<std::remove_cvref_t<T>>; ++i) {
            err = std::max(ce_abs(1.0 - approx[i] / actual[i]), err);
        }
        return err;
    } else {
        return ce_abs(1.0 - approx / actual);
    }
}

template<typename T> PF_C20CONSTEXPR double relative_l2_norm(const T &approx, const T &actual) {
    // Helper: squared absolute value, constexpr-safe for real and complex.
    const auto ce_norm = [](const auto &v) constexpr noexcept -> double {
        if constexpr (is_complex_v<std::remove_cvref_t<decltype(v)>>) {
            return double(v.real()) * double(v.real()) + double(v.imag()) * double(v.imag());
        } else {
            return double(v) * double(v);
        }
    };
    double num = 0.0;
    double denom = 0.0;
    if constexpr (has_tuple_size_v<T>) {
        for (std::size_t i = 0; i < std::tuple_size_v<std::remove_cvref_t<T>>; ++i) {
            num += ce_norm(approx[i] - actual[i]);
            denom += ce_norm(actual[i]);
        }
    } else {
        num += ce_norm(approx - actual);
        denom += ce_norm(actual);
    }
    const double ratio = denom == 0.0 ? num : num / denom;
    if (PF_IS_CONSTANT_EVALUATED()) return ce_sqrt(ratio);
    return std::sqrt(ratio);
}
} // namespace poly_eval::detail

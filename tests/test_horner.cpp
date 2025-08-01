// test_nd_horner_typed.cpp
// -----------------------------------------------------------------------------
// Random-coefficient multivariate polynomial evaluation tests for both float
// and double. Compares poly_eval::horner (run-time & compile-time) against a
// Vandermonde baseline, and exercises horner_transposed in scalar and SIMD
// modes, with both runtime and compile-time sizes.
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

#include "polyfit/fast_eval.hpp"

using std::size_t;

// templated tolerance ~100× machine epsilon
template <typename T> constexpr T eps = std::numeric_limits<T>::epsilon() * T(100);

// RNG helper
static std::mt19937 rng(42);
static std::uniform_real_distribution<double> uni_dist(-1.0, 1.0);

template <typename T> T uni() { return static_cast<T>(uni_dist(rng)); }

template <typename T> std::vector<T> random_vector(size_t n) {
    std::vector<T> v(n);
    for (auto &x : v)
        x = uni<T>();
    return v;
}

// build all multi-indices (n0,...,n_{Dim-1}), 0 ≤ ni < Deg
template <size_t Dim> auto multi_indices(size_t Deg) {
    size_t M = 1;
    for (size_t i = 0; i < Dim; ++i)
        M *= Deg;
    std::vector<std::array<size_t, Dim>> out(M);
    for (size_t k = 0; k < M; ++k) {
        size_t t = k;
        for (size_t d = 0; d < Dim; ++d) {
            out[k][d] = t % Deg;
            t /= Deg;
        }
    }
    return out;
}

// Vandermonde baseline with reversed-degree layout
template <typename T, size_t Dim>
std::array<T, Dim> vander_eval(const std::array<T, Dim> &x, const std::vector<T> &coeffs, size_t Deg) {
    auto monoms = multi_indices<Dim>(Deg);
    T powers[Dim][33]{};
    for (size_t d = 0; d < Dim; ++d) {
        powers[d][0] = T(1);
        for (size_t k = 1; k < Deg; ++k)
            powers[d][k] = powers[d][k - 1] * x[d];
    }
    auto coeff_off = [Deg](auto const &n, size_t od) {
        size_t off = od, stride = Dim;
        for (size_t d = Dim; d-- > 0;) {
            size_t rev = Deg - 1 - n[d];
            off += rev * stride;
            stride *= Deg;
        }
        return off;
    };
    std::array<T, Dim> res{};
    for (auto &n : monoms) {
        auto mono = T(1);
        for (size_t d = 0; d < Dim; ++d)
            mono *= powers[d][n[d]];
        for (size_t od = 0; od < Dim; ++od)
            res[od] += mono * coeffs[coeff_off(n, od)];
    }
    return res;
}

// random coefficients [Deg]^Dim × Dim
template <typename T, size_t Dim> auto random_coeffs(size_t Deg) {
    size_t M = 1;
    for (size_t i = 0; i < Dim; ++i)
        M *= Deg;
    M *= Dim;
    return random_vector<T>(M);
}

// unified ND Horner tester
template <typename T, size_t Dim, size_t Deg> void run_nd_horner() {
    auto coeffs = random_coeffs<T, Dim>(Deg);
    std::array<size_t, Dim + 1> exts;
    for (size_t d = 0; d < Dim; ++d)
        exts[d] = Deg;
    exts[Dim] = Dim;
    auto C = stdex::mdspan<const T, stdex::dextents<size_t, Dim + 1>>{coeffs.data(), exts};

    std::array<T, Dim> x;
    for (auto &xi : x)
        xi = uni<T>();

    auto rt = poly_eval::horner<0, true, std::array<T, Dim>>(x, C, Deg);
    auto ct = poly_eval::horner<Deg, true, std::array<T, Dim>>(x, C, Deg);
    auto vd = vander_eval<T, Dim>(x, coeffs, Deg);
    for (size_t i = 0; i < Dim; ++i) {
        EXPECT_NEAR(rt[i], vd[i], eps<T>);
        EXPECT_NEAR(rt[i], ct[i], eps<T>);
    }
}

// simple scalar Horner reference
template <typename T> T naive_horner_scalar(T x, const T *c, size_t n) {
    T acc = c[0];
    for (size_t i = 1; i < n; ++i)
        acc = acc * x + c[i];
    return acc;
}

// element-wise reference for horner_transposed
template <typename T> T naive_horner_transposed_elem(const T *x, const T *c, size_t M, size_t N, size_t i) {
    size_t stride = M;
    T acc = c[0 * stride + i];
    for (size_t k = 1; k < N; ++k)
        acc = acc * x[i] + c[k * stride + i];
    return acc;
}

// type-parameterized test suite for float/double
template <typename T> class HornerTyped : public testing::Test {};
typedef testing::Types<float, double> FloatingTypes;
TYPED_TEST_SUITE(HornerTyped, FloatingTypes);

// ND Horner sweeps
TYPED_TEST(HornerTyped, ND_Dim1_Deg2to32) {
    using T = TypeParam;
    poly_eval::detail::unroll_loop<2, 32>([](auto D) { run_nd_horner<T, 1, decltype(D)::value>(); });
}
TYPED_TEST(HornerTyped, ND_Dim2_Deg2to16) {
    using T = TypeParam;
    poly_eval::detail::unroll_loop<2, 16>([](auto D) { run_nd_horner<T, 2, decltype(D)::value>(); });
}
TYPED_TEST(HornerTyped, ND_Dim3_Deg2to8) {
    using T = TypeParam;
    poly_eval::detail::unroll_loop<2, 8>([](auto D) { run_nd_horner<T, 3, decltype(D)::value>(); });
}
TYPED_TEST(HornerTyped, ND_Dim4_Deg2to4) {
    using T = TypeParam;
    poly_eval::detail::unroll_loop<2, 4>([](auto D) { run_nd_horner<T, 4, decltype(D)::value>(); });
}

// Scalar Horner
TYPED_TEST(HornerTyped, ScalarHorner_Runtime) {
    using T = TypeParam;
    std::vector<T> c = random_vector<T>(5);
    T x = uni<T>();
    T rt = poly_eval::horner(x, c.data(), c.size());
    T ex = naive_horner_scalar(x, c.data(), c.size());
    EXPECT_NEAR(rt, ex, eps<T>);
}
TYPED_TEST(HornerTyped, ScalarHorner_CompileTime) {
    using T = TypeParam;
    constexpr size_t N = 5;
    std::vector<T> c = random_vector<T>(N);
    T x = uni<T>();
    T ct = poly_eval::horner<N>(x, c.data());
    T ex = naive_horner_scalar(x, c.data(), N);
    EXPECT_NEAR(ct, ex, eps<T>);
}

// SIMD Horner runtime
TYPED_TEST(HornerTyped, SIMDHorner_Runtime) {
    using T = TypeParam;
    std::vector<T> pts = random_vector<T>(6);
    std::vector<T> c = random_vector<T>(3);
    std::vector<T> out(pts.size());
    poly_eval::horner<0, false, false, 0>(pts.data(), out.data(), pts.size(), c.data(), c.size(),
                                          [](auto v) { return v; });
    for (size_t i = 0; i < pts.size(); ++i) {
        T ex = naive_horner_scalar(pts[i], c.data(), c.size());
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}

// horner_many runtime & compile-time
TYPED_TEST(HornerTyped, HornerMany_Runtime) {
    using T = TypeParam;
    T x = uni<T>();
    constexpr size_t M = 2, N = 3;
    std::vector<T> c = random_vector<T>(M * N);
    T out[M];
    poly_eval::horner_many<0, 0, false, T, T>(x, c.data(), out, M, N);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_scalar(x, c.data() + i * N, N);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}
TYPED_TEST(HornerTyped, HornerMany_CompileTime) {
    using T = TypeParam;
    T x = uni<T>();
    constexpr size_t M = 2, N = 2;
    std::vector<T> vec = random_vector<T>(M * N);
    T c[M * N];
    std::copy(vec.begin(), vec.end(), c);
    T out[M];
    poly_eval::horner_many<M, N, false, T, T>(x, c, out);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_scalar(x, c + i * N, N);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}

// horner_transposed scalar & SIMD, runtime & compile-time
TYPED_TEST(HornerTyped, HornerTransposed_Scalar_Runtime) {
    using T = TypeParam;
    constexpr size_t M = 6, N = 4;
    std::array<T, M> x{};
    for (auto &xi : x)
        xi = uni<T>();
    std::vector<T> c = random_vector<T>(M * N);
    std::array<T, M> out{};
    poly_eval::horner_transposed<0, 0, 0, T, T>(x.data(), c.data(), out.data(), M, N);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_transposed_elem(x.data(), c.data(), M, N, i);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}
TYPED_TEST(HornerTyped, HornerTransposed_Scalar_CompileTime) {
    using T = TypeParam;
    constexpr size_t M = 6, N = 4;
    std::array<T, M> x{};
    for (auto &xi : x)
        xi = uni<T>();
    std::vector<T> c = random_vector<T>(M * N);
    std::array<T, M> out{};
    poly_eval::horner_transposed<M, N, 0, T, T>(x.data(), c.data(), out.data(), 0, 0);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_transposed_elem(x.data(), c.data(), M, N, i);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}
TYPED_TEST(HornerTyped, HornerTransposed_SIMD_Runtime) {
    using T = TypeParam;
    using B = xsimd::batch<T>;
    constexpr size_t simd_w = B::size, M = simd_w * 2, N = 5;
    std::array<T, M> x{};
    for (auto &xi : x)
        xi = uni<T>();
    std::vector<T> c = random_vector<T>(M * N);
    std::array<T, M> out{};
    poly_eval::horner_transposed<0, 0, simd_w, T, T>(x.data(), c.data(), out.data(), M, N);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_transposed_elem(x.data(), c.data(), M, N, i);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}
TYPED_TEST(HornerTyped, HornerTransposed_SIMD_CompileTime) {
    using T = TypeParam;
    using B = xsimd::batch<T>;
    constexpr size_t simd_w = B::size, M = simd_w * 2, N = 5;
    std::array<T, M> x{};
    for (auto &xi : x)
        xi = uni<T>();
    std::vector<T> c = random_vector<T>(M * N);
    std::array<T, M> out{};
    poly_eval::horner_transposed<M, N, simd_w, T, T>(x.data(), c.data(), out.data(), 0, 0);
    for (size_t i = 0; i < M; ++i) {
        T ex = naive_horner_transposed_elem(x.data(), c.data(), M, N, i);
        EXPECT_NEAR(out[i], ex, eps<T>);
    }
}

// main test main
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

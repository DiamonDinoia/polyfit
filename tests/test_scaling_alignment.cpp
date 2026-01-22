#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "polyfit/fast_eval.hpp"

using std::size_t;

template <typename T> static T naive_horner_scalar(T x, const T *c, size_t n) {
    T acc = c[0];
    for (size_t i = 1; i < n; ++i)
        acc = acc * x + c[i];
    return acc;
}

TEST(HornerMany, ScalingTruePerPoly) {
    using T = double;
    constexpr size_t M = 2, N = 3; // two polynomials, quadratic

    // Define coefficients in reversed order per polynomial: [c2, c1, c0]
    // p0(t) = 1 + 2 t + 3 t^2
    // p1(t) = -1 + 0.5 t
    T coeffs[M * N] = {
        3.0, 2.0, 1.0, // poly 0
        0.0, 0.5, -1.0 // poly 1
    };

    // Per-poly domain mapping [a_i, b_i]
    T a[M] = {-1.0, 2.0};
    T b[M] = {+1.0, 6.0};
    T low[M], hi[M];
    for (size_t i = 0; i < M; ++i) {
        low[i] = T(1) / (b[i] - a[i]);
        hi[i] = b[i] + a[i];
    }

    // Evaluate at some x in both domains
    T x = 0.3; // in [-1,1] for poly0; maps to t in [-1,1] for poly1 via its own scaling
    T out[M]{};
    poly_eval::horner_many<M, N, true, T, T>(x, coeffs, out, 0, 0, low, hi);

    // Expected: compute each with its scaled t_i
    T t0 = (2 * x - hi[0]) * low[0];
    T t1 = (2 * x - hi[1]) * low[1];
    T ex0 = naive_horner_scalar(t0, coeffs + 0 * N, N);
    T ex1 = naive_horner_scalar(t1, coeffs + 1 * N, N);
    EXPECT_NEAR(out[0], ex0, 1e-12);
    EXPECT_NEAR(out[1], ex1, 1e-12);
}

TEST(FuncEval, AlignmentVarianceBulkEval) {
    using T = double;
    auto f = [](T x) { return std::sin(x) + T(0.25) * x; };
    T a = -1.0, b = 1.0;
    auto fe = poly_eval::make_func_eval(f, 8, a, b);

    const size_t P = 65; // use odd + offsets to vary alignment
    std::vector<T> pts(P + 2), out1(P + 2), out2(P + 2);
    for (size_t i = 0; i < P + 2; ++i)
        pts[i] = a + (b - a) * (T(i) / T(P + 1));

    // Mismatched alignments (different offsets) to exercise safe unaligned path
    fe(pts.data() + 1, out2.data() + 0, P);
    for (size_t i = 0; i < P; ++i) {
        EXPECT_NEAR(out2[i], fe(pts[i + 1]), 1e-12);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

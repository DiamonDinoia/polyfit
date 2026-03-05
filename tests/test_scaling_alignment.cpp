#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "polyfit/fast_eval.hpp"

using std::size_t;

template<typename T> static T naive_horner_scalar(T x, const T *c, size_t n) {
    T acc = c[0];
    for (size_t i = 1; i < n; ++i) acc = acc * x + c[i];
    return acc;
}

TEST(HornerMany, ScalingPerPoly) {
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

    // Evaluate at some x — apply per-poly domain mapping, then call horner per poly
    T x = 0.3;
    T out[M]{};
    for (size_t i = 0; i < M; ++i) {
        T xm = (T(2) * x - hi[i]) * low[i];
        out[i] = naive_horner_scalar(xm, coeffs + i * N, N);
    }

    // Verify against horner_many (no scaling) with pre-mapped x per poly
    for (size_t i = 0; i < M; ++i) {
        T xm = (T(2) * x - hi[i]) * low[i];
        T result[1]{};
        poly_eval::horner_many<1, N, T, T>(xm, coeffs + i * N, result);
        EXPECT_NEAR(result[0], out[i], 1e-12);
    }
}

TEST(FuncEval, AlignmentVarianceBulkEval) {
    using T = double;
    auto f = [](T x) {
        return std::sin(x) + 0.25 * x;
    };
    T a = -1.0, b = 1.0;
    auto fe = poly_eval::make_func_eval(f, 8, a, b);

    const size_t P = 65; // use odd + offsets to vary alignment
    std::vector<T> pts(P + 2), out1(P + 2), out2(P + 2);
    for (size_t i = 0; i < P + 2; ++i) pts[i] = a + (b - a) * (T(i) / T(P + 1));

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

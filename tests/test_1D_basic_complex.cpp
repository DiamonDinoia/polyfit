// test_1D_basic_complex.cpp: runtime/compile-time complex<double> fitting tests

#include "test_1D_helpers.hpp"

static std::mt19937 gen(42);

// 5. Runtime Coefficient Count with complex<double>
TEST(PolyEval, RuntimeDegreeComplexRandom) {
    double a = -1.0, b = 1.0;
    int nCoeffs = 13;
    auto poly = poly_eval::fit(complex_func, nCoeffs, a, b);
    const auto eps = 1e-6;
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}

TEST(PolyEval, CompileDegreeComplexRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr auto N = 13;
    const auto poly = poly_eval::fit<N>(complex_func, a, b);
    const auto eps = 1e-12;
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, eps);
}

TEST(PolyEval, RuntimeEpsComplexRandom) {
    double a = -1.0, b = 1.0;
    double eps = 1e-13;
    auto poly = poly_eval::fit(complex_func, eps, a, b);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}

// test_1D_truncation.cpp — truncation, adaptive fit, and coefficient-order tests

#include "test_1D_helpers.h"

static std::mt19937 gen(42);

TEST(PolyEval, TruncateLowDegreeFunc) {
    auto cubic = [](double x) { return x * x * x - 2.0 * x + 1.0; };
    auto poly = poly_eval::fit(cubic, 16, -1.0, 1.0);
    const auto original_size = poly.coeffs().size();
    EXPECT_EQ(original_size, 16u);

    poly.truncate(1e-8);
    EXPECT_LT(poly.coeffs().size(), original_size);
    EXPECT_LE(poly.coeffs().size(), 6u);

    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t i = 0; i < 100; ++i) {
        double x = dist(gen);
        EXPECT_NEAR(poly(x), cubic(x), 1e-7) << "x=" << x;
    }
}

TEST(PolyEval, TruncatePreservesConstant) {
    auto constantFunc = [](double x) { (void)x; return 42.0; };
    auto poly = poly_eval::fit(constantFunc, 16, -1.0, 1.0);
    poly.truncate(1e-10);
    EXPECT_EQ(poly.coeffs().size(), 1u);
    EXPECT_NEAR(poly(0.5), 42.0, 1e-10);
}

TEST(PolyEval, CoeffsAreHighDegreeFirst) {
    auto quadratic = [](double x) { return x * x + 2.0; };
    auto poly = poly_eval::fit(quadratic, 3, -1.0, 1.0);
    const auto &coeffs = poly.coeffs();
    ASSERT_EQ(coeffs.size(), 3u);
    EXPECT_NEAR(coeffs[0], 1.0, 1e-12);
    EXPECT_NEAR(coeffs[1], 0.0, 1e-12);
    EXPECT_NEAR(coeffs[2], 2.0, 1e-12);
}

TEST(PolyEval, AdaptiveFitFindsWorkingDegree) {
    auto poly = poly_eval::fit(double_func, 1e-12, -0.5, 0.5);
    EXPECT_LE(poly.coeffs().size(), 32u);
    EXPECT_GT(poly.coeffs().size(), 1u);

    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (std::size_t i = 0; i < 100; ++i) {
        double x = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(poly(x), double_func(x)), 1e-12);
    }
}

// Edge-case and special floating-point value tests for polyfit
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

#include "polyfit/fast_eval.hpp"

// --------------------------------------------------------------------------
// Edge cases: batch evaluation with num_points=0 and num_points=1
// --------------------------------------------------------------------------

static auto cos_func = [](double x) { return std::cos(x); };

TEST(EdgeCases, BatchEvalZeroPoints) {
    auto poly = poly_eval::make_func_eval(cos_func, 8, -1.0, 1.0);
    std::vector<double> xs;
    std::vector<double> ys;
    // Should be a no-op without crashing
    poly(xs.data(), ys.data(), 0);
}

TEST(EdgeCases, BatchEvalSinglePoint) {
    auto poly = poly_eval::make_func_eval(cos_func, 16, -1.0, 1.0);
    double x = 0.5;
    double y = 0.0;
    poly(&x, &y, 1);
    EXPECT_NEAR(y, std::cos(0.5), 1e-12);
}

// --------------------------------------------------------------------------
// Edge case: domain where a == b
// --------------------------------------------------------------------------

TEST(EdgeCases, DegenerateDomainAEqualsB) {
    // When a == b the domain mapping involves division by zero (1/(b-a)).
    // The evaluator constructs but produces Inf/NaN on evaluation.
    auto identity = [](double x) { return x; };
    auto poly = poly_eval::make_func_eval(identity, 4, 1.0, 1.0);
    const double result = poly(1.0);
    EXPECT_TRUE(std::isinf(result) || std::isnan(result))
        << "Expected inf or NaN for degenerate domain, got " << result;
}

// --------------------------------------------------------------------------
// Edge case: very small domain
// --------------------------------------------------------------------------

TEST(EdgeCases, VerySmallDomain) {
    constexpr double center = 1.0;
    constexpr double half = 1e-12;
    auto poly = poly_eval::make_func_eval(cos_func, 8, center - half, center + half);
    const double mid = center;
    EXPECT_NEAR(poly(mid), std::cos(mid), 1e-14);
}

// --------------------------------------------------------------------------
// NaN / Infinity behavior documentation tests
// --------------------------------------------------------------------------

TEST(SpecialValues, NaNInputProducesNaN) {
    auto poly = poly_eval::make_func_eval(cos_func, 8, -1.0, 1.0);
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    const double result = poly(nan_val);
    // NaN propagation: evaluating at NaN should produce NaN
    EXPECT_TRUE(std::isnan(result)) << "Expected NaN, got " << result;
}

TEST(SpecialValues, InfinityInputProducesInfOrNaN) {
    auto poly = poly_eval::make_func_eval(cos_func, 8, -1.0, 1.0);
    const double inf_val = std::numeric_limits<double>::infinity();
    const double result = poly(inf_val);
    // Evaluating at infinity should produce inf or NaN (not a finite wrong answer)
    EXPECT_TRUE(std::isinf(result) || std::isnan(result))
        << "Expected inf or NaN, got " << result;
}

TEST(SpecialValues, NaNBatchInputProducesNaN) {
    auto poly = poly_eval::make_func_eval(cos_func, 8, -1.0, 1.0);
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    double out = 0.0;
    poly(&nan_val, &out, 1);
    EXPECT_TRUE(std::isnan(out)) << "Expected NaN in batch output, got " << out;
}

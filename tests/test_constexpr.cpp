// test_constexpr.cpp
// Static-assert coverage for compile-time evaluation paths:
//   detail::math::{fma,abs,log10,sqrt} and FuncEval CT construction.
// All key tests are static_assert — no runtime state needed.

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

#include "polyfit/polyfit.hpp"
#include "polyfit/internal/macros.h"
#include "polyfit/internal/utils.h"

namespace m = poly_eval::detail::math;

#if PF_HAS_CXX20 || (!defined(_MSC_VER) && (defined(__clang__) || defined(__GNUC__)))
#define PF_TEST_CONSTEXPR_MATH 1
#else
#define PF_TEST_CONSTEXPR_MATH 0
#endif

// ---------------------------------------------------------------------------
// math::fma
// ---------------------------------------------------------------------------
#if PF_TEST_CONSTEXPR_MATH
static_assert(m::fma(2.0, 3.0, 1.0) == 7.0, "fma(2,3,1)");
static_assert(m::fma(0.0, 1.0, 5.0) == 5.0, "fma(0,1,5)");
static_assert(m::fma(-1.0, 2.0, 4.0) == 2.0, "fma(-1,2,4)");
static_assert(m::fma(2.0f, 3.0f, 1.0f) == 7.0f, "fma float");

// ---------------------------------------------------------------------------
// math::abs
// ---------------------------------------------------------------------------
static_assert(m::abs(-3.14) > 3.0, "abs(-3.14) > 3");
static_assert(m::abs(2.5) == 2.5, "abs(2.5)");
static_assert(m::abs(-2.5) == 2.5, "abs(-2.5)");
static_assert(m::abs(0.0) == 0.0, "abs(0)");
static_assert(m::abs(-1.0f) == 1.0f, "abs float");

// ---------------------------------------------------------------------------
// math::log10
// ---------------------------------------------------------------------------
static_assert(m::log10(1.0) == 0.0, "log10(1) == 0");
static_assert(m::log10(10.0) > 0.99 && m::log10(10.0) < 1.01, "log10(10) ~= 1");
static_assert(m::log10(100.0) > 1.99 && m::log10(100.0) < 2.01, "log10(100) ~= 2");
static_assert(m::log10(0.1) > -1.01 && m::log10(0.1) < -0.99, "log10(0.1) ~= -1");

// ---------------------------------------------------------------------------
// math::sqrt
// ---------------------------------------------------------------------------
static_assert(m::sqrt(4.0) == 2.0, "sqrt(4)");
static_assert(m::sqrt(0.0) == 0.0, "sqrt(0)");
static_assert(m::sqrt(9.0) == 3.0, "sqrt(9)");
static_assert(m::sqrt(2.0) > 1.41 && m::sqrt(2.0) < 1.42, "sqrt(2) ~= 1.414"); // NOLINT(modernize-use-std-numbers)
#endif

// ---------------------------------------------------------------------------
// FuncEval CT construction (C++20+)
// ---------------------------------------------------------------------------
#if __cplusplus >= 202002L
static constexpr auto ct_poly = poly_eval::fit<4>(
    [](double x) { return x * x + 1.0; }, -1.0, 1.0);

static_assert(decltype(ct_poly)::NCOEFFS == 4, "FuncEval CT: NCOEFFS");
static_assert(poly_eval::detail::relativeL2Norm(ct_poly(0.0), 1.0) < 1e-12,
              "FuncEval CT: eval at x=0 should approximate x^2+1=1");

static constexpr auto ct_fused_poly = poly_eval::fit<3>(
    [](double x) constexpr { return x * x + 2.0 * x + 1.0; }, 2.0, 4.0, poly_eval::FuseAlways{});
static constexpr auto ct_fused_coeffs = ct_fused_poly.coeffs();

static_assert(poly_eval::detail::relativeL2Norm(ct_fused_poly(3.0), 16.0) < 1e-12,
              "FuncEval CT FuseAlways: evaluator should use fused coefficients on the original domain");
static_assert(poly_eval::detail::relativeL2Norm(
                  poly_eval::horner(3.0, ct_fused_coeffs.data(), ct_fused_coeffs.size()), 16.0) < 1e-12,
              "FuncEval CT FuseAlways: coeffs() should match runtime Horner evaluation on the original domain");

static constexpr std::array<double, 2> nd_a{-1.0, -1.0};
static constexpr std::array<double, 2> nd_b{1.0, 1.0};
static constexpr auto ct_nd = poly_eval::fit<4, nd_a, nd_b>([](const std::array<double, 2> &p) constexpr {
    return std::array<double, 2>{p[0] + p[1], p[0] - p[1]};
});

static_assert(ct_nd.nCoeffsPerAxis() == 4, "FuncEvalND CT: nCoeffsPerAxis");
#endif

// ---------------------------------------------------------------------------
// Minimal runtime smoke test
// ---------------------------------------------------------------------------
TEST(Constexpr, StaticAssertsPassed) {
    SUCCEED();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#undef PF_TEST_CONSTEXPR_MATH

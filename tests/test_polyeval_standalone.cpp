// Verifies that <polyfit/polyeval.hpp> is usable on its own. The
// evaluation-only API must not pull in fitting machinery, and `horner` /
// `hybrid` must produce the documented results on user-supplied coefficients.

#include <polyfit/polyeval.hpp>

#include <gtest/gtest.h>

#include <array>

namespace {

// p(x) = 2*x^3 - 3*x^2 + 4*x - 5, coefficients in Horner order
// (highest-degree first).
constexpr std::array<double, 4> kCoeffs{2.0, -3.0, 4.0, -5.0};

constexpr double referencePoly(double x) noexcept {
    return ((2.0 * x - 3.0) * x + 4.0) * x - 5.0;
}

TEST(PolyevalStandalone, HornerCompileTimeDegreeMatchesReference) {
    for (double x : {-2.0, -0.5, 0.0, 0.25, 1.0, 3.5}) {
        const double got = poly_eval::horner<kCoeffs.size()>(x, kCoeffs.data());
        EXPECT_DOUBLE_EQ(got, referencePoly(x));
    }
}

TEST(PolyevalStandalone, HornerRuntimeDegreeMatchesReference) {
    for (double x : {-2.0, -0.5, 0.0, 0.25, 1.0, 3.5}) {
        const double got = poly_eval::horner<>(x, kCoeffs.data(), kCoeffs.size());
        EXPECT_DOUBLE_EQ(got, referencePoly(x));
    }
}

TEST(PolyevalStandalone, HybridMatchesHorner) {
    for (double x : {-2.0, -0.5, 0.0, 0.25, 1.0, 3.5}) {
        const double horner = poly_eval::horner<kCoeffs.size()>(x, kCoeffs.data());
        const double hybrid = poly_eval::hybrid<kCoeffs.size()>(x, kCoeffs.data());
        EXPECT_DOUBLE_EQ(hybrid, horner);
    }
}

} // namespace

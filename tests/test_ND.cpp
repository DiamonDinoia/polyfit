#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <random>

#include "polyfit/polyfit.hpp"

// -----------------------------------------------------------------------------
// Test settings
// -----------------------------------------------------------------------------

constexpr int kNumPoints = 1000; // random evaluation points per test

// Global test function: sum of cosines across all dimensions.  Replace or edit
// the body if you want to test a different analytic function.
template<typename Array, typename OutPut = Array> constexpr auto sumCos(const Array &x) {
    double s = 0.0;
    for (const double xi : x) s += std::cos(xi);
    OutPut y{};
    y.fill(s);
    for (std::size_t i = 1; i < y.size(); ++i) {
        y[i] += y[i - 1];
    }
    return y;
}

std::mt19937 gen(42);
std::uniform_real_distribution<double> dist(-1.0, 1.0);

// -----------------------------------------------------------------------------
// Helper template that runs a complete Chebyshev correctness test for the
// compile‑time combination <InDim,OutDim,nCoeffs>.
// -----------------------------------------------------------------------------

template<std::size_t inDim, std::size_t outDim, std::size_t nCoeffsCt> void runMonomialTest() {
    using Input = std::array<double, inDim>;
    using Output = std::array<double, outDim>;
    const double tol = std::pow(10.0, -static_cast<double>(nCoeffsCt - 3));
    Input a{};
    a.fill(-1.0);
    Input b{};
    b.fill(1.0);

    auto approx = poly_eval::make_func_eval<nCoeffsCt>(sumCos<Input, Output>, a, b);

    for (int i = 0; i < kNumPoints; ++i) {
        Input x;
        for (auto &xi : x) xi = dist(gen);
        auto expected = sumCos<Input, Output>(x);
        auto actual = approx(x);
        for (std::size_t j = 0; j < outDim; ++j) {
            ASSERT_NEAR(actual[j], expected[j], tol);
        }
    }
}

// -----------------------------------------------------------------------------
// Individual TEST cases – no custom preprocessor tricks, just standard GTest.
// -----------------------------------------------------------------------------

TEST(Eval, In2Out2Deg16) { runMonomialTest<2, 2, 16>(); }
TEST(Eval, In2Out3Deg16) { runMonomialTest<2, 3, 16>(); }
TEST(Eval, In3Out2Deg16) { runMonomialTest<3, 2, 16>(); }
TEST(Eval, In3Out3Deg8) { runMonomialTest<3, 3, 8>(); }
TEST(Eval, In3Out4Deg8) { runMonomialTest<3, 4, 8>(); }
TEST(Eval, In4Out3Deg8) { runMonomialTest<4, 3, 8>(); }
TEST(Eval, In4Out4Deg4) { runMonomialTest<4, 4, 4>(); }

TEST(Eval, RuntimeNCoeffsRejectsNonPositive) {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;
    In a{-1.0, -1.0};
    In b{1.0, 1.0};
    EXPECT_THROW((void)poly_eval::make_func_eval(sumCos<In, Out>, 0, a, b), std::invalid_argument);
    EXPECT_THROW((void)poly_eval::make_func_eval(sumCos<In, Out>, -2, a, b), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Main entry point for Google Test
// -----------------------------------------------------------------------------
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

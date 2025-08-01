#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <random>

#include "polyfit/fast_eval.hpp"

// -----------------------------------------------------------------------------
// Test settings
// -----------------------------------------------------------------------------

constexpr int kNumPoints = 1000; // random evaluation points per test

// Global test function: sum of cosines across all dimensions.  Replace or edit
// the body if you want to test a different analytic function.
template <typename Array, typename OutPut = Array> constexpr auto sumCos(Array const &x) {
    double s = 0.0;
    for (const double xi : x)
        s += std::cos(xi);
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
// compile‑time combination <InDim,OutDim,Degree>.
// -----------------------------------------------------------------------------

template <std::size_t InDim, std::size_t OutDim, std::size_t Degree> void RunmMonomialTest() {
    using Input = std::array<double, InDim>;
    using Output = std::array<double, OutDim>;
    // set ktol to 10^{-Degree}
    constexpr double kTol = std::pow(10.0, -static_cast<double>(Degree - 3));
    // Domain [a,b] = [-1,1]^InDim
    Input a{};
    a.fill(-1.0);
    Input b{};
    b.fill(1.0);

    // Build Chebyshev approximation once
    auto approx = poly_eval::make_func_eval<Degree>(sumCos<Input, Output>, a, b);

    // Pseudo‑random validation points
    for (int i = 0; i < kNumPoints; ++i) {
        Input x;
        for (auto &xi : x)
            xi = dist(gen);
        auto y_true = sumCos<Input, Output>(x);
        auto y_poly = approx(x);
        for (std::size_t j = 0; j < OutDim; ++j) {
            ASSERT_NEAR(y_poly[j], y_true[j], kTol);
        }
    }
}

// -----------------------------------------------------------------------------
// Individual TEST cases – no custom preprocessor tricks, just standard GTest.
// -----------------------------------------------------------------------------

TEST(Eval, In2_Out2_Deg16) { RunmMonomialTest<2, 2, 16>(); }
TEST(Eval, In2_Out3_Deg16) { RunmMonomialTest<2, 3, 16>(); }
TEST(Eval, In3_Out2_Deg16) { RunmMonomialTest<3, 2, 16>(); }
TEST(Eval, In3_Out3_Deg8) { RunmMonomialTest<3, 3, 8>(); }
TEST(Eval, In3_Out4_Deg8) { RunmMonomialTest<3, 4, 8>(); }
TEST(Eval, In4_Out3_Deg8) { RunmMonomialTest<4, 3, 8>(); }
TEST(Eval, In4_Out4_Deg2) { RunmMonomialTest<4, 4, 4>(); }

// -----------------------------------------------------------------------------
// Main entry point for Google Test
// -----------------------------------------------------------------------------
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

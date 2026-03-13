#include <gtest/gtest.h>
#include <random>

#include "test_ND_shared.hpp"

std::mt19937 gen(42);
std::uniform_real_distribution<double> dist(-1.0, 1.0);

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS> void runMonomialTest() {
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;
    const double tol = std::pow(10.0, -static_cast<double>(NCOEFFS - 3));
    Input a{};
    a.fill(-1.0);
    Input b{};
    b.fill(1.0);

    auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, a, b);

    for (int i = 0; i < kNumPoints; ++i) {
        Input x;
        for (auto &xi : x) xi = dist(gen);
        auto expected = sumCos<Input, Output>(x);
        auto actual = approx(x);
        for (std::size_t j = 0; j < OUT_DIM; ++j) ASSERT_NEAR(actual[j], expected[j], tol);
    }
}

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
    const In a{-1.0, -1.0};
    const In b{1.0, 1.0};
    EXPECT_THROW((void)poly_eval::fit(sumCos<In, Out>, 0, a, b), std::invalid_argument);
    EXPECT_THROW((void)poly_eval::fit(sumCos<In, Out>, -2, a, b), std::invalid_argument);
}

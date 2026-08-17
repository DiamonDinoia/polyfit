// Monomial accuracy TESTs live one-instantiation-per-TU in
// test_ND_accuracy_in*_out*_deg*.cpp; this TU keeps the cheap API checks.

#include <gtest/gtest.h>

#include "test_ND_shared.hpp"

TEST(Eval, RuntimeNCoeffsRejectsNonPositive) {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;
    const In a{-1.0, -1.0};
    const In b{1.0, 1.0};
    EXPECT_THROW((void)poly_eval::fit(sumCos<In, Out>, 0, a, b), std::invalid_argument);
    EXPECT_THROW((void)poly_eval::fit(sumCos<In, Out>, -2, a, b), std::invalid_argument);
}

// polyeval.hpp undefines the PF_* macros at its end; polyfit.hpp included afterwards must
// still see them, so the macro headers cannot be include-guarded.

#include <polyfit/polyeval.hpp>
#include <polyfit/polyfit.hpp>

#include <gtest/gtest.h>

#include <array>

TEST(IncludeOrder, PolyevalThenPolyfit) {
    constexpr std::array<double, 3> c{1.0, -2.0, 1.0}; // (x - 1)^2
    EXPECT_EQ(poly_eval::horner<3>(3.0, c.data()), 4.0);
    const auto fit = poly_eval::fit<4>([](double x) { return x * x; }, -1.0, 1.0);
    EXPECT_NEAR(fit(0.5), 0.25, 1e-14);
}

// test_ND_fusion_wide_domain.cpp — FusionMode::Always accuracy off-centre.

#include <algorithm>
#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

#include "test_ND_fusion_shared.hpp"

using namespace polyfit_test_nd_fusion;

TEST(PolyEvalND, FusionAlwaysOnWideDomainStillMeetsTol) {
    // Off-centre, non-unit domain — fusion must fold alpha/beta in and the
    // forced Always path must keep within a modest tolerance of the exact
    // function on random interior points.
    auto f = [](const Arr2 &x) -> Out1 {
        return {std::exp(-x[0]) * std::cos(x[1])};
    };
    const Arr2 a{-0.25, 0.5};
    const Arr2 b{1.75, 2.25};

    auto poly_force =
        poly_eval::FuncEvalND<decltype(f), 14, poly_eval::FusionMode::Always>(f, a, b);

    std::mt19937 gen(7);
    std::uniform_real_distribution<double> dx(a[0], b[0]);
    std::uniform_real_distribution<double> dy(a[1], b[1]);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        Arr2 x{dx(gen), dy(gen)};
        const double approx = poly_force(x)[0];
        const double exact = f(x)[0];
        maxErr = std::max(maxErr, std::abs(approx - exact));
    }
    EXPECT_LT(maxErr, 1e-10);
}

// test_ND_fusion_matches_auto.cpp — FusionMode::Always vs Auto agreement.

#include <random>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

#include "test_ND_fusion_shared.hpp"

using namespace polyfit_test_nd_fusion;

TEST(PolyEvalND, FusionAlwaysMatchesAutoOnNarrowDomain) {
    auto poly_always =
        poly_eval::FuncEvalND<decltype(smooth_2d), 10, poly_eval::FusionMode::Always>(
            smooth_2d, Arr2{-1.0, -1.0}, Arr2{1.0, 1.0});
    auto poly_auto =
        poly_eval::FuncEvalND<decltype(smooth_2d), 10, poly_eval::FusionMode::Auto>(
            smooth_2d, Arr2{-1.0, -1.0}, Arr2{1.0, 1.0});

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        Arr2 x{d(gen), d(gen)};
        const double ya = poly_always(x)[0];
        const double yb = poly_auto(x)[0];
        EXPECT_NEAR(ya, yb, 1e-10);
    }
}

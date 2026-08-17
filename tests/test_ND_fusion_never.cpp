// test_ND_fusion_never.cpp — FusionMode::Never skips domain folding.

#include <algorithm>
#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

#include "test_ND_fusion_shared.hpp"

using namespace polyfit_test_nd_fusion;

TEST(PolyEvalND, FusionNeverSkipsFusion) {
    auto poly_never =
        poly_eval::FuncEvalND<decltype(smooth_2d), 10, poly_eval::FusionMode::Never>(
            smooth_2d, Arr2{-1.0, -1.0}, Arr2{1.0, 1.0});

    std::mt19937 gen(5);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        Arr2 x{d(gen), d(gen)};
        maxErr = std::max(maxErr, std::abs(poly_never(x)[0] - smooth_2d(x)[0]));
    }
    EXPECT_LT(maxErr, 1e-8);
}

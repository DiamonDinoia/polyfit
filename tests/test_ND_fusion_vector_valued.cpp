// test_ND_fusion_vector_valued.cpp: FusionMode::Always with multi-output f.

#include <algorithm>
#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

#include "test_ND_fusion_shared.hpp"

using namespace polyfit_test_nd_fusion;

TEST(PolyEvalND, FusionAlwaysVectorValued) {
    auto f = [](const Arr2 &x) -> Out2 {
        return {std::sin(x[0] + x[1]), std::cos(x[0] - x[1])};
    };

    auto poly =
        poly_eval::FuncEvalND<decltype(f), 12, poly_eval::FusionMode::Always>(
            f, Arr2{0.0, 0.0}, Arr2{1.0, 1.0});

    std::mt19937 gen(11);
    std::uniform_real_distribution<double> d(0.0, 1.0);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        Arr2 x{d(gen), d(gen)};
        const auto approx = poly(x);
        const auto exact = f(x);
        for (std::size_t k = 0; k < 2; ++k)
            maxErr = std::max(maxErr, std::abs(approx[k] - exact[k]));
    }
    EXPECT_LT(maxErr, 1e-10);
}

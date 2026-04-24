// test_ND_fusion.cpp — FuncEvalND fusion (per-axis) tests.

#include <array>
#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <polyfit/polyeval.hpp>

namespace {
constexpr std::size_t kNumRandomTests = 4000;

using Arr2 = std::array<double, 2>;
using Out1 = std::array<double, 1>;
using Out2 = std::array<double, 2>;

auto smooth_2d = [](const Arr2 &x) -> Out1 {
    return {std::sin(1.3 * x[0]) * std::cos(0.7 * x[1])};
};

} // namespace

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

TEST(PolyEvalND, FusionAlwaysDoesNotStoreDomainParams) {
    using FE_fused =
        poly_eval::FuncEvalND<decltype(smooth_2d), 8, poly_eval::FusionMode::Always>;
    using FE_never =
        poly_eval::FuncEvalND<decltype(smooth_2d), 8, poly_eval::FusionMode::Never>;
    // With FusionMode::Always the per-axis domain params and identity flag
    // should collapse into the [[no_unique_address]] empty storage.
    static_assert(sizeof(FE_fused) < sizeof(FE_never),
                  "Always evaluator should be smaller than Never (no domain params stored)");
}

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

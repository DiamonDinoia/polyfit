// test_1D_fusion.cpp — FuseNever/FuseAlways/tag-order tests

#include "test_1D_helpers.h"

static std::mt19937 gen(42);

TEST(PolyEval, FusionNeverSkipsFusion) {
    auto func = [](double x) { return std::sin(x); };
    const double a = 0.0, b = poly_eval::detail::constants::pi;
    auto poly_never = poly_eval::fit<16>(func, a, b, poly_eval::FuseNever{});
    auto poly_auto  = poly_eval::fit<16>(func, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr_never = 0.0, maxErr_auto = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        maxErr_never = std::max(maxErr_never, poly_eval::detail::relativeL2Norm(poly_never(x), func(x)));
        maxErr_auto  = std::max(maxErr_auto,  poly_eval::detail::relativeL2Norm(poly_auto(x),  func(x)));
    }
    EXPECT_LT(maxErr_never, 1e-10);
    EXPECT_LT(maxErr_auto,  1e-10);
    EXPECT_NE(maxErr_never, maxErr_auto)
        << "FuseNever and auto should take different code paths on non-trivial domains";
}

TEST(PolyEval, FusionAlwaysForcesOnNarrowDomain) {
    auto func = [](double x) { return std::cos(x); };
    auto poly_always = poly_eval::fit<8>(func, 0.0, 2.0, poly_eval::FuseAlways{});
    auto poly_auto   = poly_eval::fit<8>(func, 0.0, 2.0);

    std::uniform_real_distribution<double> dist(0.0, 2.0);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        EXPECT_DOUBLE_EQ(poly_always(x), poly_auto(x)) << "x=" << x;
    }
}

TEST(PolyEval, FusionAlwaysOnWideDomainLosesAccuracy) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    auto poly_auto  = poly_eval::fit(func, 32, a, b);
    auto poly_force = poly_eval::fit(func, 32, a, b, poly_eval::FuseAlways{});

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr_auto = 0.0, maxErr_force = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        maxErr_auto  = std::max(maxErr_auto,  poly_eval::detail::relativeL2Norm(poly_auto(x),  func(x)));
        maxErr_force = std::max(maxErr_force, poly_eval::detail::relativeL2Norm(poly_force(x), func(x)));
    }
    EXPECT_LT(maxErr_auto, 1e-13);
    EXPECT_GT(maxErr_force, maxErr_auto)
        << "auto: " << maxErr_auto << ", force: " << maxErr_force;
}

TEST(PolyEval, FuseAlwaysDoesNotStoreDomainParams) {
    auto func = [](double x) { return x * x; };
    using FE_fused = decltype(poly_eval::fit<4>(func, -1.0, 1.0, poly_eval::FuseAlways{}));
    using FE_never = decltype(poly_eval::fit<4>(func, -1.0, 1.0, poly_eval::FuseNever{}));
    // FuseAlways should not store invSpan, sumEndpoints, identityMap
    static_assert(sizeof(FE_fused) < sizeof(FE_never),
                  "FuseAlways evaluator should be smaller than FuseNever (no domain params stored)");
}

TEST(PolyEval, TagOrderIndependence) {
    auto func = [](double x) { return std::sin(x); };
    auto p1 = poly_eval::fit<16>(func, -1.0, 1.0, poly_eval::Iters<2>{}, poly_eval::FuseNever{});
    auto p2 = poly_eval::fit<16>(func, -1.0, 1.0, poly_eval::FuseNever{}, poly_eval::Iters<2>{});

    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        EXPECT_DOUBLE_EQ(p1(x), p2(x)) << "x=" << x;
    }
}

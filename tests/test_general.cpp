#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <random>

#include "polyfit/polyfit.hpp"

std::mt19937 rng(42);

static_assert(poly_eval::detail::FitOptions<>::VALID);
static_assert(poly_eval::detail::FitOptions<poly_eval::EvalPts<64>, poly_eval::MaxCoeffs<12>, poly_eval::Iters<2>,
                                            poly_eval::FuseNever>::VALID);
static_assert(poly_eval::detail::FitOptions<poly_eval::MaxCoeffs<8>, poly_eval::MaxCoeffs<16>>::VALID == false);
static_assert(poly_eval::detail::FitOptions<poly_eval::FuseAuto, poly_eval::FuseNever>::VALID == false);

// Maximum relative error over random samples of a 2D domain.
template<typename TrueF, typename ApproxF, typename Domain>
double computeMaxRelativeError(TrueF f_true, ApproxF f_approx, Domain low, Domain high, int num_samples) {
    std::uniform_real_distribution<double> dist_x(low[0], high[0]);
    std::uniform_real_distribution<double> dist_y(low[1], high[1]);

    double maxErr = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        const Domain pt{dist_x(rng), dist_y(rng)};
        double err = poly_eval::detail::relativeError(f_approx(pt), f_true(pt));
        maxErr = std::max(maxErr, err);
    }
    return maxErr;
}

TEST(FuncEval2D, FixedNCoeffsVsAdaptiveRelativeError) {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    auto f = [](const In &in) {
        auto [x, y] = in;
        return Out{std::cos(3 * x + y) * std::sin(x - 2 * y), std::exp(1 + x * x + y * y)};
    };
    In low{-1.0, -1.0};
    In high{1.0, 1.0};
    constexpr int num_samples = 1000;

    // Fixed coefficient count
    constexpr int nCoeffsFixed = 16;
    auto fe_fixed = poly_eval::fit(f, nCoeffsFixed, low, high);
    // Adaptive
    double eps = 1e-7;
    auto fe_adapt = poly_eval::fit(f, eps, low, high);

    double maxRelFixed = computeMaxRelativeError(
        [&](const In &pt) { return f(pt); }, [&](const In &pt) { return fe_fixed(pt); }, low, high, num_samples);
    double maxRelAdapt = computeMaxRelativeError(
        [&](const In &pt) { return f(pt); }, [&](const In &pt) { return fe_adapt(pt); }, low, high, num_samples);

    EXPECT_LT(maxRelFixed, 1e-5);
    EXPECT_LT(maxRelAdapt, 1e-4);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

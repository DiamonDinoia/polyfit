#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <tuple>

#include "polyfit/fast_eval.hpp"

std::mt19937 rng(42);

// Helper to compute maximum relative error over random samples
template <typename TrueF, typename ApproxF, typename Domain>
double compute_max_relative_error(TrueF f_true, ApproxF f_approx, Domain low, Domain high, int num_samples) {
    // RNG seeded from GoogleTest's random seed for reproducibility

    std::uniform_real_distribution<double> dist_x;
    std::uniform_real_distribution<double> dist_y;

    if constexpr (std::is_floating_point_v<Domain>) {
        dist_x = std::uniform_real_distribution<double>(low, high);
    } else {
        dist_x = std::uniform_real_distribution<double>(std::get<0>(low), std::get<0>(high));
        dist_y = std::uniform_real_distribution<double>(std::get<1>(low), std::get<1>(high));
    }

    double max_err = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        Domain pt{};
        if constexpr (std::is_floating_point_v<Domain>) {
            pt = dist_x(rng);
        } else {
            pt = Domain{dist_x(rng), dist_y(rng)};
        }
        const auto true_val = f_true(pt);
        const auto approx_val = f_approx(pt);
        double err = poly_eval::detail::relative_error(approx_val, true_val);
        max_err = std::max(max_err, err);
    }
    return max_err;
}

TEST(FuncEval1D, FixedDegreeVsAdaptiveRelativeError) {
    using namespace poly_eval;
    auto f = [](double x) { return std::exp(std::sin(3 * x)); };
    double a = -1.0;
    double b = 1.0;
    constexpr int num_samples = 1000;

    // Fixed-degree
    constexpr int N_fixed = 16;
    auto fe_fixed = poly_eval::make_func_eval(f, N_fixed, a, b);
    // Adaptive
    constexpr double eps = 1e-6;
    auto fe_adapt = poly_eval::make_func_eval(f, eps, a, b);

    double max_rel_fixed = compute_max_relative_error([&](double x) { return f(x); },
                                                      [&](double x) { return fe_fixed(x); }, a, b, num_samples);
    double max_rel_adapt = compute_max_relative_error([&](double x) { return f(x); },
                                                      [&](double x) { return fe_adapt(x); }, a, b, num_samples);

    EXPECT_LT(max_rel_fixed, 1e-2);
    EXPECT_LT(max_rel_adapt, 1e-5);
}

TEST(FuncEval2D, FixedDegreeVsAdaptiveRelativeError) {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    auto f = [](const In &in) {
        auto [x, y] = in;
        return Out{std::cos(3 * x + y) * std::sin(x - 2 * y), std::exp(1 + x * x + y * y)};
    };
    In low{-1.0, -1.0};
    In high{1.0, 1.0};
    constexpr int num_samples = 1000;

    // Fixed-degree
    constexpr int N_fixed = 16;
    auto fe_fixed = poly_eval::make_func_eval(f, N_fixed, low, high);
    // Adaptive
    double eps = 1e-7;
    auto fe_adapt = poly_eval::make_func_eval(f, eps, low, high);

    double max_rel_fixed = compute_max_relative_error(
        [&](const In &pt) { return f(pt); }, [&](const In &pt) { return fe_fixed(pt); }, low, high, num_samples);
    double max_rel_adapt = compute_max_relative_error(
        [&](const In &pt) { return f(pt); }, [&](const In &pt) { return fe_adapt(pt); }, low, high, num_samples);

    EXPECT_LT(max_rel_fixed, 1e-5);
    EXPECT_LT(max_rel_adapt, 1e-4);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

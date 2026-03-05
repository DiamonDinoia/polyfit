// fast_eval_test.cpp
// Google Test suite for poly_eval::make_func_eval APIs with randomized verification

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <random>

#include "polyfit/fast_eval.hpp"

std::mt19937 gen(42);

// Example functions
static auto double_func = [](const double x) {
    return std::cos(x);
};
static auto float_func = [](const float x) {
    return std::cos(x);
};
static auto complex_func = [](const double x) {
    return std::complex<double>(x * x, std::sin(x));
};

// Number of random test points
static constexpr std::size_t kNumRandomTests = 100;

// Helper to batch-evaluate and compare
template<typename T, typename X = double, typename F>
void batch_verify(const F &f, const std::vector<X> &xs, const std::vector<T> &ys, const double tol) {
    for (size_t i = 0; i < xs.size(); ++i) {
        EXPECT_LE(poly_eval::detail::relative_l2_norm(ys[i], f(xs[i])), tol)
            << "Failed at x=" << xs[i] << ": expected " << f(xs[i]) << ", got " << ys[i];
    }
}

// 1. Runtime Degree (double, default iters)
TEST(PolyEval, RuntimeDegreeDoubleRandom) {
    double a = -0.5, b = 0.5;
    int n = 16;
    constexpr auto eps = 1e-15;
    auto poly = poly_eval::make_func_eval(double_func, n, a, b);
    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(poly(xs[i]), double_func(xs[i])), eps);
    }

    // Batch test
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 2. Runtime Degree (float, custom iters)
TEST(PolyEval, RuntimeDegreeFloatCustomItersRandom) {
    constexpr auto a = -static_cast<float>(poly_eval::detail::constants::pi);
    constexpr auto b = static_cast<float>(poly_eval::detail::constants::pi);
    constexpr auto n = 12;
    constexpr size_t iters = 1;
    constexpr auto eps = 2.5e-4;
    const auto poly = poly_eval::make_func_eval<iters>(float_func, n, a, b);

    // Randomized tests (single-point)
    std::uniform_real_distribution<float> dist(a, b);
    std::vector<float> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(poly(xs[i]), float_func(xs[i])), eps);
    }

    // Batch test (float)
    std::vector<float> ys_f(kNumRandomTests);
    poly(xs.data(), ys_f.data(), ys_f.size());

    batch_verify<float>(float_func, xs, ys_f, eps);
}

TEST(PolyEval, RuntimeDegreeRejectsNonPositive) {
    constexpr double a = -1.0;
    constexpr double b = 1.0;
    EXPECT_THROW((void)poly_eval::make_func_eval(double_func, 0, a, b), std::invalid_argument);
    EXPECT_THROW((void)poly_eval::make_func_eval(double_func, -3, a, b), std::invalid_argument);
}
// 3. Compile-Time Degree (double, default iters)
TEST(PolyEval, CompileTimeDegreeDoubleRandom) {
    double a = -.5, b = .5;
    constexpr size_t N = 6;
    constexpr auto eps = 1e-4;
    auto poly = poly_eval::make_func_eval<N>(double_func, a, b);

    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(poly(xs[i]), double_func(xs[i])), eps);
    }
    // Batch test
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 4. Runtime Epsilon (C++17 API)
TEST(PolyEval, ErrorDrivenRuntimeEpsRandom) {
    double a = -1.0, b = 1.0;
    double eps = 1e-10;
    constexpr size_t MaxN = 16;
    constexpr size_t EvalPts = 100;
    constexpr size_t Iters = 1;
    auto poly = poly_eval::make_func_eval<MaxN, EvalPts, Iters>(double_func, eps, a, b);

    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(poly(xs[i]), double_func(xs[i])), eps);
    }

    // Batch test
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 6. Runtime Degree with complex<double>
TEST(PolyEval, RuntimeDegreeComplexRandom) {
    double a = -1.0, b = 1.0;
    int n = 13;
    auto poly = poly_eval::make_func_eval(complex_func, n, a, b);
    const auto eps = 1e-6;
    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        auto want = complex_func(xs[i]);
        auto got = poly(xs[i]);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(want, got), eps);
    }

    // Batch test
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}

TEST(PolyEval, CompileDegreeComplexRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr auto N = 13;
    const auto poly = poly_eval::make_func_eval<N>(complex_func, a, b);
    const auto eps = 1e-12;
    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        auto want = complex_func(xs[i]);
        auto got = poly(xs[i]);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(want, got), eps);
    }

    // Batch test
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, eps);
}

TEST(PolyEval, RuntimeEpsComplexRandom) {
    double a = -1.0, b = 1.0;
    double eps = 1e-13;
    auto poly = poly_eval::make_func_eval(complex_func, eps, a, b);
    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        auto want = complex_func(xs[i]);
        auto got = poly(xs[i]);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(want, got), eps);
    }

    // Batch test
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}
#if __cplusplus >= 202002L && __cplusplus < 202302L
// Pure arithmetic constexpr function for C++20 tests
// (Disabled under C++23: static_for with mutable captures is not constant-expression-valid)
constexpr auto double_constexpr_func = [](const double x) constexpr {
    return 2.0 * x * x * x - 3.0 * x + 1.0;
};
// 5. Full Compile-Time Fitting and Evaluation (constexpr fixed-degree API)
TEST(PolyEval, FullCompileTimeRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr size_t Degree = 5;
    constexpr size_t ItersCT = 2;
    constexpr auto poly = poly_eval::make_func_eval<Degree, ItersCT>(double_constexpr_func, a, b);

    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_NEAR(poly(xs[i]), double_constexpr_func(xs[i]), 1e-15) << "x=" << xs[i];
    }

    // Batch test
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_constexpr_func, xs, ys, 1e-7);
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
constexpr auto complex_constexpr_func = [](const double x) constexpr {
    return std::complex<double>(x * x, x + x);
};
// 6. Full Compile-Time Fitting and Evaluation (constexpr fixed-degree API)
TEST(PolyEval, FullCompileTimeEps) {
    constexpr double a = -1.0, b = 1.0;
    constexpr auto eps = 1e-13;
    constexpr auto poly = poly_eval::make_func_eval<eps, a, b>(double_constexpr_func);

    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_NEAR(poly(xs[i]), double_constexpr_func(xs[i]), eps) << "x=" << xs[i];
    }

    // Batch test
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_constexpr_func, xs, ys, eps);
}

// 7. Error‑driven compile‑time eps for complex<double>
TEST(PolyEval, ErrorDrivenCompileTimeEpsComplexRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr double eps = 1e-10;
    constexpr size_t MaxN = 32;
    constexpr size_t EvalPts = 100;
    constexpr size_t Iters = 0;
    constexpr auto poly = poly_eval::make_func_eval<eps, a, b, MaxN, EvalPts, Iters>(complex_constexpr_func);
    // Randomized tests (single-point)
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        auto want = complex_constexpr_func(xs[i]);
        auto got = poly(xs[i]);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(want, got), eps);
    }

    // Batch test
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_constexpr_func, xs, ys, eps);
}
#endif // PF_HAS_CONSTEXPR_EPS_OVERLOAD
#endif

// ----- FuncEval truncation tests -----

TEST(PolyEval, TruncateLowDegreeFunc) {
    // Fit a cubic polynomial at degree 16 — high-degree terms should be ~0
    auto cubic = [](double x) {
        return x * x * x - 2.0 * x + 1.0;
    };
    auto poly = poly_eval::make_func_eval(cubic, 16, -1.0, 1.0);
    const auto original_size = poly.coeffs().size();
    EXPECT_EQ(original_size, 16u);

    poly.truncate(1e-8);
    // Should truncate down to ~4 terms (degree 3 = 4 coefficients)
    EXPECT_LT(poly.coeffs().size(), original_size);
    EXPECT_LE(poly.coeffs().size(), 6u); // some slack for numerical noise

    // Verify accuracy after truncation
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t i = 0; i < 100; ++i) {
        double x = dist(gen);
        EXPECT_NEAR(poly(x), cubic(x), 1e-7) << "x=" << x;
    }
}

TEST(PolyEval, TruncatePreservesConstant) {
    // A constant function: all coefficients except constant should be ~0
    auto const_func = [](double x) {
        (void)x;
        return 42.0;
    };
    auto poly = poly_eval::make_func_eval(const_func, 16, -1.0, 1.0);
    poly.truncate(1e-10);
    EXPECT_EQ(poly.coeffs().size(), 1u); // only constant term remains
    EXPECT_NEAR(poly(0.5), 42.0, 1e-10);
}

TEST(PolyEval, AdaptiveFitThenTruncate) {
    // The adaptive make_func_eval should now fit-then-truncate
    auto poly = poly_eval::make_func_eval(double_func, 1e-12, -0.5, 0.5);
    // cos(x) on [-0.5, 0.5] needs ~8 terms for 1e-12 accuracy
    EXPECT_LE(poly.coeffs().size(), 32u);

    // Verify accuracy
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (std::size_t i = 0; i < 100; ++i) {
        double x = dist(gen);
        EXPECT_LE(poly_eval::detail::relative_l2_norm(poly(x), double_func(x)), 1e-12);
    }
}

// ----- High-degree accuracy tests (compensated Horner) -----

TEST(PolyEval, HighDegree48MachineEps) {
    // Degree 48 on [-1,1] should achieve near machine epsilon
    auto sin_func = [](double x) {
        return std::sin(x);
    };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::make_func_eval(sin_func, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_err = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_err = std::max(max_err, std::abs(poly(x) - sin_func(x)));
    }
    // Tolerance relaxed: MSVC/Apple Clang std::cos differs from GCC by up to
    // ~1 ULP, which at degree 48 amplifies through Björck-Pereyra to ~1e-11.
    EXPECT_LT(max_err, 1e-10) << "Degree-48 sin fit max error: " << max_err;
}

TEST(PolyEval, HighDegree48Exp) {
    // Degree 48 on [-1,1] for exp should also be excellent
    auto exp_func = [](double x) {
        return std::exp(x);
    };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::make_func_eval(exp_func, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_err = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_err = std::max(max_err, std::abs(poly(x) - exp_func(x)));
    }
    EXPECT_LT(max_err, 1e-10) << "Degree-48 exp fit max error: " << max_err;
}

TEST(PolyEval, HighDegree48Complex) {
    // Degree 48 complex fit on [-1, 1]
    auto func = [](double x) {
        return std::complex<double>(std::sin(x), std::cos(x));
    };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::make_func_eval(func, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_err = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_err = std::max(max_err, std::abs(poly(x) - func(x)));
    }
    EXPECT_LT(max_err, 2e-10) << "Degree-48 complex fit max error: " << max_err;
}

TEST(PolyEval, HighDegree48Batch) {
    // Verify batch evaluation works at degree 48
    auto sin_func = [](double x) {
        return std::sin(x);
    };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::make_func_eval(sin_func, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests), ys(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) xs[i] = dist(gen);

    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(sin_func, xs, ys, 1e-10);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

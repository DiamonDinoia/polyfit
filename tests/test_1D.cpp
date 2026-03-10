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
    const auto poly = poly_eval::make_func_eval(float_func, n, a, b, poly_eval::iters<iters>{});

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
    auto poly = poly_eval::make_func_eval(double_func, eps, a, b, poly_eval::max_degree<MaxN>{},
                                          poly_eval::eval_pts<EvalPts>{}, poly_eval::iters<Iters>{});

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
    constexpr auto poly = poly_eval::make_func_eval<Degree>(double_constexpr_func, a, b, poly_eval::iters<ItersCT>{});

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
    constexpr auto poly =
        poly_eval::make_func_eval<eps, a, b, MaxN, EvalPts, Iters>(complex_constexpr_func);
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

// ----- Wide-offset domain tests (domain mapping accuracy) -----
// Fitting on [1000, 2000] exercises the affine map to [-1, 1]:
//   x_mapped = (2*x - (a+b)) / (b-a)
// Catastrophic cancellation in the map can lose significant digits when
// |a+b| >> |b-a|, so these tests verify the full pipeline preserves accuracy.
// Functions must be smooth enough on the interval for moderate-degree approximation.

TEST(PolyEval, WideDomainRuntimeDegree) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::make_func_eval(func, 32, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relative_l2_norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain [1000,2000] exp(-x/1000) max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainCompileTimeDegree) {
    auto func = [](double x) { return std::log(x); };
    constexpr double a = 1000.0, b = 2000.0;
    constexpr std::size_t N = 24;
    auto poly = poly_eval::make_func_eval<N>(func, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relative_l2_norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain CT log(x) max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainBatch) {
    auto func = [](double x) { return 1.0 / x; };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::make_func_eval(func, 32, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests), ys(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) xs[i] = dist(gen);

    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(func, xs, ys, 1e-13);
}

TEST(PolyEval, WideDomainEpsDriven) {
    auto func = [](double x) { return std::sqrt(x); };
    const double a = 1000.0, b = 2000.0;
    constexpr double eps = 1e-12;
    auto poly = poly_eval::make_func_eval(func, eps, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relative_l2_norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, eps) << "Wide domain eps-driven sqrt max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainComplex) {
    auto func = [](double x) {
        return std::complex<double>(std::log(x), 1.0 / x);
    };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::make_func_eval(func, 32, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relative_l2_norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain complex log+i/x max relative error: " << max_rel;
}

// ----- Direct-domain interpolation (no [-1,1] mapping) comparison -----
// Demonstrates the accuracy catastrophe when interpolating directly on a
// wide-offset domain like [1000, 2000] without mapping to [-1, 1].
// The monomial basis {1, x, x^2, ...} with x ~ 1500 has condition number
// ~ 1500^deg, destroying all precision at moderate degrees.

TEST(PolyEval, DirectDomainVsMappedAccuracy) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    constexpr std::size_t N = 16;

    // --- Approach 1: with domain mapping (the normal API) ---
    auto poly_mapped = poly_eval::make_func_eval(func, static_cast<int>(N), a, b);

    // --- Approach 2: direct interpolation on [a, b] without mapping ---
    // Place Chebyshev nodes directly in [a, b]
    std::vector<double> nodes(N), samples(N);
    for (std::size_t k = 0; k < N; ++k) {
        double t = std::cos((2.0 * double(k) + 1.0) * poly_eval::detail::constants::pi / (2.0 * double(N)));
        nodes[k] = 0.5 * ((b - a) * t + (b + a)); // Chebyshev node in [a, b]
        samples[k] = func(nodes[k]);
    }
    auto newton = poly_eval::detail::bjorck_pereyra(nodes, samples);
    auto mono = poly_eval::detail::newton_to_monomial(newton, nodes);

    // Evaluate both at random points
    std::uniform_real_distribution<double> dist(a, b);
    double max_err_mapped = 0.0, max_err_direct = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        double exact = func(x);
        double err_mapped = std::abs(poly_mapped(x) - exact);
        double err_direct = std::abs(poly_eval::horner(x, mono.data(), mono.size()) - exact);
        max_err_mapped = std::max(max_err_mapped, err_mapped);
        max_err_direct = std::max(max_err_direct, err_direct);
    }

    // The mapped version should be accurate
    EXPECT_LT(max_err_mapped, 1e-13)
        << "Mapped max abs error: " << max_err_mapped;

    // The direct (unmapped) version should be MUCH worse — many orders of magnitude
    // At degree 16 with nodes at ~1500, monomial conditioning ~ 1500^16 ~ 10^50,
    // so we expect total loss of precision.
    EXPECT_GT(max_err_direct, max_err_mapped * 1e3)
        << "Direct interpolation should be much worse than mapped.\n"
        << "  mapped max error:  " << max_err_mapped << "\n"
        << "  direct max error:  " << max_err_direct;

    // Print for human inspection
    std::printf("\n  [DirectDomainVsMapped] N=%zu, domain=[%.0f, %.0f]\n", N, a, b);
    std::printf("    Mapped ([-1,1]) max abs error: %.2e\n", max_err_mapped);
    std::printf("    Direct (no map) max abs error: %.2e\n", max_err_direct);
    std::printf("    Ratio direct/mapped: %.1fx worse\n", max_err_direct / max_err_mapped);
}

// ----- Fusion mode tests -----

TEST(PolyEval, FusionNeverSkipsFusion) {
    // Use [0, pi] where the domain map is non-trivial (not identity like [-1,1])
    auto func = [](double x) { return std::sin(x); };
    const double a = 0.0, b = poly_eval::detail::constants::pi;
    auto poly_never = poly_eval::make_func_eval<16>(func, a, b, poly_eval::fuse_never{});
    auto poly_auto = poly_eval::make_func_eval<16>(func, a, b);
    std::uniform_real_distribution<double> dist(a, b);
    double max_err_never = 0.0, max_err_auto = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_err_never = std::max(max_err_never, poly_eval::detail::relative_l2_norm(poly_never(x), func(x)));
        max_err_auto = std::max(max_err_auto, poly_eval::detail::relative_l2_norm(poly_auto(x), func(x)));
    }
    EXPECT_LT(max_err_never, 1e-10);
    EXPECT_LT(max_err_auto, 1e-10);
    EXPECT_NE(max_err_never, max_err_auto)
        << "fuse_never and auto should take different code paths on non-trivial domains";
}

TEST(PolyEval, FusionAlwaysForcesOnNarrowDomain) {
    auto func = [](double x) { return std::cos(x); };
    auto poly_always = poly_eval::make_func_eval<8>(func, 0.0, 2.0, poly_eval::fuse_always{});
    auto poly_auto = poly_eval::make_func_eval<8>(func, 0.0, 2.0);
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        EXPECT_DOUBLE_EQ(poly_always(x), poly_auto(x)) << "x=" << x;
    }
}

TEST(PolyEval, FusionAlwaysOnWideDomainLosesAccuracy) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    auto poly_auto = poly_eval::make_func_eval(func, 32, a, b);
    auto poly_force = poly_eval::make_func_eval(func, 32, a, b, poly_eval::fuse_always{});

    std::uniform_real_distribution<double> dist(a, b);
    double max_err_auto = 0.0, max_err_force = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_err_auto = std::max(max_err_auto, poly_eval::detail::relative_l2_norm(poly_auto(x), func(x)));
        max_err_force = std::max(max_err_force, poly_eval::detail::relative_l2_norm(poly_force(x), func(x)));
    }
    EXPECT_LT(max_err_auto, 1e-13);
    EXPECT_GT(max_err_force, max_err_auto)
        << "auto: " << max_err_auto << ", force: " << max_err_force;
}

TEST(PolyEval, TagOrderIndependence) {
    auto func = [](double x) { return std::sin(x); };
    auto p1 = poly_eval::make_func_eval<16>(func, -1.0, 1.0, poly_eval::iters<2>{}, poly_eval::fuse_never{});
    auto p2 = poly_eval::make_func_eval<16>(func, -1.0, 1.0, poly_eval::fuse_never{}, poly_eval::iters<2>{});
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        EXPECT_DOUBLE_EQ(p1(x), p2(x)) << "x=" << x;
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

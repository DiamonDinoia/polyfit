// test_1D_basic.cpp — runtime/compile-time degree and eps-driven fitting tests

#include "test_1D_helpers.h"

static std::mt19937 gen(42);

struct EvalCarrier {
    double value = 0.0;

    constexpr EvalCarrier() = default;
    constexpr EvalCarrier(double v) : value(v) {}
};

constexpr EvalCarrier operator*(EvalCarrier lhs, EvalCarrier rhs) noexcept { return EvalCarrier(lhs.value * rhs.value); }
constexpr EvalCarrier operator+(EvalCarrier lhs, EvalCarrier rhs) noexcept { return EvalCarrier(lhs.value + rhs.value); }

// 1. Runtime Coefficient Count (double, default iters)
TEST(PolyEval, RuntimeDegreeDoubleRandom) {
    double a = -0.5, b = 0.5;
    int nCoeffs = 16;
    constexpr auto eps = 1e-15;
    auto poly = poly_eval::fit(double_func, nCoeffs, a, b);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(poly(xs[i]), double_func(xs[i])), eps);
    }
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 2. Runtime Coefficient Count (float, custom iters)
TEST(PolyEval, RuntimeDegreeFloatCustomItersRandom) {
    constexpr auto a = -static_cast<float>(poly_eval::detail::constants::pi);
    constexpr auto b =  static_cast<float>(poly_eval::detail::constants::pi);
    constexpr auto nCoeffs = 12;
    constexpr size_t iters = 1;
    constexpr auto eps = 2.5e-4;
    const auto poly = poly_eval::fit(float_func, nCoeffs, a, b, poly_eval::Iters<iters>{});
    std::uniform_real_distribution<float> dist(a, b);
    std::vector<float> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(poly(xs[i]), float_func(xs[i])), eps);
    }
    std::vector<float> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), ys.size());
    batch_verify<float>(float_func, xs, ys, eps);
}

TEST(PolyEval, RuntimeDegreeRejectsNonPositive) {
    constexpr double a = -1.0;
    constexpr double b =  1.0;
    EXPECT_THROW((void)poly_eval::fit(double_func, 0,  a, b), std::invalid_argument);
    EXPECT_THROW((void)poly_eval::fit(double_func, -3, a, b), std::invalid_argument);
}

// 3. Compile-Time Coefficient Count (double, default iters)
TEST(PolyEval, CompileTimeDegreeDoubleRandom) {
    double a = -.5, b = .5;
    constexpr size_t N = 6;
    constexpr auto eps = 1e-4;
    auto poly = poly_eval::fit<N>(double_func, a, b);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(poly(xs[i]), double_func(xs[i])), eps);
    }
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 4. Runtime Epsilon (C++17 API)
TEST(PolyEval, ErrorDrivenRuntimeEpsRandom) {
    double a = -1.0, b = 1.0;
    double eps = 1e-10;
    constexpr size_t maxNCoeffs = 16;
    constexpr size_t evalPoints = 100;
    constexpr size_t Iters = 1;
    auto poly = poly_eval::fit(double_func, eps, a, b, poly_eval::MaxCoeffs<maxNCoeffs>{},
                               poly_eval::EvalPts<evalPoints>{}, poly_eval::Iters<Iters>{});
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(poly(xs[i]), double_func(xs[i])), eps);
    }
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_func, xs, ys, eps);
}

// 5. Runtime Coefficient Count with complex<double>
TEST(PolyEval, RuntimeDegreeComplexRandom) {
    double a = -1.0, b = 1.0;
    int nCoeffs = 13;
    auto poly = poly_eval::fit(complex_func, nCoeffs, a, b);
    const auto eps = 1e-6;
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}

TEST(PolyEval, CompileDegreeComplexRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr auto N = 13;
    const auto poly = poly_eval::fit<N>(complex_func, a, b);
    const auto eps = 1e-12;
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, eps);
}

TEST(PolyEval, RuntimeEpsComplexRandom) {
    double a = -1.0, b = 1.0;
    double eps = 1e-13;
    auto poly = poly_eval::fit(complex_func, eps, a, b);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_func, xs, ys, 1e-6);
}

#if __cplusplus >= 202002L && __cplusplus < 202302L
constexpr auto double_constexpr_func = [](const double x) constexpr {
    return 2.0 * x * x * x - 3.0 * x + 1.0;
};

TEST(PolyEval, FullCompileTimeRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr size_t nCoeffs = 5;
    constexpr size_t ItersCT = 2;
    constexpr auto poly = poly_eval::fit<nCoeffs>(double_constexpr_func, a, b, poly_eval::Iters<ItersCT>{});
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_NEAR(poly(xs[i]), double_constexpr_func(xs[i]), 1e-15) << "x=" << xs[i];
    }
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_constexpr_func, xs, ys, 1e-7);
}

TEST(PolyEval, GenericScalarOperatorAcceptsBatchCarrier) {
    constexpr double a = -1.0;
    constexpr double b = 1.0;
    constexpr std::size_t nCoeffs = 6;
    const auto poly = poly_eval::fit<nCoeffs>(double_func, a, b);

    using Batch = xsimd::batch<double>;
    alignas(Batch::arch_type::alignment()) std::array<double, Batch::size> xs{};
    alignas(Batch::arch_type::alignment()) std::array<double, Batch::size> ys{};
    for (std::size_t i = 0; i < xs.size(); ++i) xs[i] = a + (b - a) * (static_cast<double>(i) / static_cast<double>(xs.size()));

    const auto xb = Batch::load_aligned(xs.data());
    const auto yb = poly(xb);
    yb.store_aligned(ys.data());

    for (std::size_t i = 0; i < xs.size(); ++i) {
        EXPECT_NEAR(ys[i], poly(xs[i]), 1e-15);
    }
}

TEST(PolyEval, GenericScalarOperatorAcceptsCustomCarrier) {
    constexpr double a = -1.0;
    constexpr double b = 1.0;
    constexpr std::size_t nCoeffs = 6;
    const auto poly = poly_eval::fit<nCoeffs>(double_func, a, b, poly_eval::FuseAlways{});

    const EvalCarrier x(0.25);
    const auto y = poly(x);
    EXPECT_NEAR(y.value, poly(x.value), 1e-15);
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
constexpr auto complex_constexpr_func = [](const double x) constexpr {
    return std::complex<double>(x * x, x + x);
};

TEST(PolyEval, FullCompileTimeEps) {
    constexpr double a = -1.0, b = 1.0;
    constexpr auto eps = 1e-13;
    constexpr auto poly = poly_eval::fit<eps, a, b>(double_constexpr_func);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_NEAR(poly(xs[i]), double_constexpr_func(xs[i]), eps) << "x=" << xs[i];
    }
    std::vector<double> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(double_constexpr_func, xs, ys, eps);
}

TEST(PolyEval, ErrorDrivenCompileTimeEpsComplexRandom) {
    constexpr double a = -1.0, b = 1.0;
    constexpr double eps = 1e-10;
    constexpr size_t maxNCoeffs = 32;
    constexpr size_t evalPoints = 100;
    constexpr size_t Iters = 0;
    constexpr auto poly =
        poly_eval::fit<eps, a, b, maxNCoeffs, evalPoints, Iters>(complex_constexpr_func);
    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        xs[i] = dist(gen);
        EXPECT_LE(poly_eval::detail::relativeL2Norm(complex_constexpr_func(xs[i]), poly(xs[i])), eps);
    }
    std::vector<std::complex<double>> ys(kNumRandomTests);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<std::complex<double>>(complex_constexpr_func, xs, ys, eps);
}
#endif // PF_HAS_CONSTEXPR_EPS_OVERLOAD
#endif // C++20

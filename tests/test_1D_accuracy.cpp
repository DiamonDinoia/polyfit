// test_1D_accuracy.cpp — high-degree accuracy and wide-domain tests

#include "test_1D_helpers.h"

static std::mt19937 gen(42);

// ----- Large-fit accuracy tests (compensated Horner) -----

TEST(PolyEval, HighDegree48MachineEps) {
    auto sinFunc = [](double x) { return std::sin(x); };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::fit(sinFunc, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        maxErr = std::max(maxErr, std::abs(poly(x) - sinFunc(x)));
    }
    EXPECT_LT(maxErr, 1e-10) << "nCoeffs-48 sin fit max error: " << maxErr;
}

TEST(PolyEval, HighDegree48Exp) {
    auto expFunc = [](double x) { return std::exp(x); };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::fit(expFunc, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        maxErr = std::max(maxErr, std::abs(poly(x) - expFunc(x)));
    }
    EXPECT_LT(maxErr, 1e-10) << "nCoeffs-48 exp fit max error: " << maxErr;
}

TEST(PolyEval, HighDegree48Complex) {
    auto func = [](double x) { return std::complex<double>(std::sin(x), std::cos(x)); };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::fit(func, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        maxErr = std::max(maxErr, std::abs(poly(x) - func(x)));
    }
    EXPECT_LT(maxErr, 2e-10) << "nCoeffs-48 complex fit max error: " << maxErr;
}

TEST(PolyEval, HighDegree48Batch) {
    auto sinFunc = [](double x) { return std::sin(x); };
    const double a = -1.0, b = 1.0;
    auto poly = poly_eval::fit(sinFunc, 48, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    std::vector<double> xs(kNumRandomTests), ys(kNumRandomTests);
    for (std::size_t i = 0; i < kNumRandomTests; ++i) xs[i] = dist(gen);
    poly(xs.data(), ys.data(), xs.size());
    batch_verify<double>(sinFunc, xs, ys, 1e-10);
}

// ----- Wide-offset domain tests (domain mapping accuracy) -----

TEST(PolyEval, WideDomainRuntimeDegree) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::fit(func, 32, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relativeL2Norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain [1000,2000] exp(-x/1000) max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainCompileTimeDegree) {
    auto func = [](double x) { return std::log(x); };
    constexpr double a = 1000.0, b = 2000.0;
    constexpr std::size_t N = 24;
    auto poly = poly_eval::fit<N>(func, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relativeL2Norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain CT log(x) max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainBatch) {
    auto func = [](double x) { return 1.0 / x; };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::fit(func, 32, a, b);

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
    auto poly = poly_eval::fit(func, eps, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relativeL2Norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, eps) << "Wide domain eps-driven sqrt max relative error: " << max_rel;
}

TEST(PolyEval, WideDomainComplex) {
    auto func = [](double x) { return std::complex<double>(std::log(x), 1.0 / x); };
    const double a = 1000.0, b = 2000.0;
    auto poly = poly_eval::fit(func, 32, a, b);

    std::uniform_real_distribution<double> dist(a, b);
    double max_rel = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x = dist(gen);
        max_rel = std::max(max_rel, poly_eval::detail::relativeL2Norm(poly(x), func(x)));
    }
    EXPECT_LT(max_rel, 1e-13) << "Wide domain complex log+i/x max relative error: " << max_rel;
}

// ----- Direct-domain interpolation comparison -----

TEST(PolyEval, DirectDomainVsMappedAccuracy) {
    auto func = [](double x) { return std::exp(-x / 1000.0); };
    const double a = 1000.0, b = 2000.0;
    constexpr std::size_t N = 16;

    auto poly_mapped = poly_eval::fit(func, static_cast<int>(N), a, b);

    std::vector<double> nodes(N), samples(N);
    for (std::size_t k = 0; k < N; ++k) {
        double t = std::cos((2.0 * double(k) + 1.0) * poly_eval::detail::constants::pi / (2.0 * double(N)));
        nodes[k]   = 0.5 * ((b - a) * t + (b + a));
        samples[k] = func(nodes[k]);
    }
    auto newton = poly_eval::detail::bjorckPereyra(nodes, samples);
    auto mono   = poly_eval::detail::newtonToMonomial(newton, nodes);

    std::uniform_real_distribution<double> dist(a, b);
    double maxErr_mapped = 0.0, maxErr_direct = 0.0;
    for (std::size_t i = 0; i < kNumRandomTests; ++i) {
        double x     = dist(gen);
        double exact = func(x);
        maxErr_mapped = std::max(maxErr_mapped, std::abs(poly_mapped(x) - exact));
        maxErr_direct = std::max(maxErr_direct,
                                 std::abs(poly_eval::horner(x, mono.data(), mono.size()) - exact));
    }

    EXPECT_LT(maxErr_mapped, 1e-13)
        << "Mapped max abs error: " << maxErr_mapped;
    EXPECT_GT(maxErr_direct, maxErr_mapped * 1e3)
        << "Direct interpolation should be much worse than mapped.\n"
        << "  mapped max error:  " << maxErr_mapped << "\n"
        << "  direct max error:  " << maxErr_direct;

    std::printf("\n  [DirectDomainVsMapped] N=%zu, domain=[%.0f, %.0f]\n", N, a, b);
    std::printf("    Mapped ([-1,1]) max abs error: %.2e\n", maxErr_mapped);
    std::printf("    Direct (no map) max abs error: %.2e\n", maxErr_direct);
    std::printf("    Ratio direct/mapped: %.1fx worse\n", maxErr_direct / maxErr_mapped);
}

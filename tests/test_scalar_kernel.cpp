// Tests for the type-level ScalarKernel selector on FuncEval / FuncEvalND.

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>

#include "polyfit/polyfit.hpp"

namespace {

template<typename T> constexpr T tol = std::numeric_limits<T>::epsilon() * T(64);

double sinFunc(double x) { return std::sin(x); }

TEST(ScalarKernel, FuncEvalDefaultIsHybrid) {
    using HybridEval = poly_eval::FuncEval<double(*)(double), 12, 1,
                                            poly_eval::FusionMode::Never>;
    using ExplicitHybrid = poly_eval::FuncEval<double(*)(double), 12, 1,
                                                poly_eval::FusionMode::Never,
                                                poly_eval::ScalarKernel::Hybrid>;
    static_assert(std::is_same_v<HybridEval, ExplicitHybrid>,
                  "Default ScalarKernel must be Hybrid");
}

TEST(ScalarKernel, FuncEvalNDDefaultIsHybrid) {
    using HybridEvalND = poly_eval::FuncEvalND<std::array<double, 2>(*)(std::array<double, 2>),
                                                8, poly_eval::FusionMode::Never>;
    using ExplicitHybridND = poly_eval::FuncEvalND<std::array<double, 2>(*)(std::array<double, 2>),
                                                    8, poly_eval::FusionMode::Never,
                                                    poly_eval::ScalarKernel::Hybrid>;
    static_assert(std::is_same_v<HybridEvalND, ExplicitHybridND>,
                  "Default ScalarKernel must be Hybrid for FuncEvalND");
}

TEST(ScalarKernel, HornerAndHybridAgreeWithin1D) {
    using HybridEv = poly_eval::FuncEval<double(*)(double), 16, 1,
                                           poly_eval::FusionMode::Never,
                                           poly_eval::ScalarKernel::Hybrid>;
    using HornerEv = poly_eval::FuncEval<double(*)(double), 16, 1,
                                           poly_eval::FusionMode::Never,
                                           poly_eval::ScalarKernel::Horner>;

    HybridEv hyb(sinFunc, -1.0, 1.0);
    HornerEv hor(sinFunc, -1.0, 1.0);

    // Coefficients are identical; the kernel choice does not affect the fit.
    ASSERT_EQ(hyb.coeffs().size(), hor.coeffs().size());
    for (std::size_t i = 0; i < hyb.coeffs().size(); ++i) {
        ASSERT_EQ(hyb.coeffs()[i], hor.coeffs()[i]) << "i=" << i;
    }

    for (double x = -1.0; x <= 1.0; x += 0.0625) {
        const double a = hyb(x);
        const double b = hor(x);
        const double ref = std::sin(x);
        EXPECT_NEAR(b, ref, tol<double>) << "Horner @x=" << x;
        EXPECT_NEAR(a, ref, tol<double>) << "Hybrid @x=" << x;
        // Within a few ulps of each other; two schemes evaluating the same
        // coefficients differ only by accumulation order.
        EXPECT_NEAR(a, b, tol<double>) << "Kernel disagreement @x=" << x;
    }
}

std::array<double, 1> ndFunc(std::array<double, 2> v) {
    return {std::sin(v[0]) * std::cos(v[1])};
}

TEST(ScalarKernel, HornerAndHybridAgreeWithinND) {
    using HybridND = poly_eval::FuncEvalND<decltype(&ndFunc), 10,
                                            poly_eval::FusionMode::Never,
                                            poly_eval::ScalarKernel::Hybrid>;
    using HornerND = poly_eval::FuncEvalND<decltype(&ndFunc), 10,
                                            poly_eval::FusionMode::Never,
                                            poly_eval::ScalarKernel::Horner>;

    const std::array<double, 2> a{-1.0, -1.0};
    const std::array<double, 2> b{1.0, 1.0};
    HybridND hyb(&ndFunc, a, b);
    HornerND hor(&ndFunc, a, b);

    // With 10 coeffs per axis over [-1,1] the fit, not the evaluator, sets the
    // accuracy. Use a loose tolerance for function-value comparison and a tight
    // one for kernel-vs-kernel agreement.
    constexpr double fitTol = 1e-6;
    for (double x = -0.9; x <= 0.9; x += 0.3) {
        for (double y = -0.9; y <= 0.9; y += 0.3) {
            const auto vh = hyb(std::array<double, 2>{x, y});
            const auto vH = hor(std::array<double, 2>{x, y});
            const double ref = std::sin(x) * std::cos(y);
            EXPECT_NEAR(vh[0], ref, fitTol);
            EXPECT_NEAR(vH[0], ref, fitTol);
            EXPECT_NEAR(vh[0], vH[0], tol<double>);
        }
    }
}

} // namespace

#include <array>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

#include "polyfit/polyfit.hpp"

using namespace poly_eval;

namespace {

template<class FE> FE rebuild1D(const FE &fe) {
    return FE(fromCoeffs, fe.coeffs(), fe.domainInvSpan(), fe.domainSumEndpoints());
}

template<class FE> FE rebuildND(const FE &fe) { return FE(fromCoeffs, fe.coeffsMdspan(), fe.domainView()); }

template<class FE> void expect1DBitExact(const FE &fe, double a, double b) {
    const FE re = rebuild1D(fe);
    for (int i = 0; i <= 64; ++i) {
        const double x = a + (b - a) * i / 64.0;
        EXPECT_EQ(fe(x), re(x)) << "x=" << x;
    }
}

template<class FE> void expectNDBitExact(const FE &fe) {
    const FE re = rebuildND(fe);
    for (int i = 0; i <= 16; ++i)
        for (int j = 0; j <= 16; ++j) {
            const std::array<double, 2> x{-1.0 + i / 8.0, 0.5 + j / 16.0};
            EXPECT_EQ(fe(x), re(x)) << "x=" << x[0] << "," << x[1];
        }
}

} // namespace

TEST(FromCoeffs, OneDimCompileTimeCount) {
    expect1DBitExact(fit<9>([](double x) { return std::exp(x); }, 0.3, 2.1), 0.3, 2.1);
}

TEST(FromCoeffs, OneDimRuntimeCount) {
    expect1DBitExact(fit([](double x) { return std::sin(3 * x); }, 11, -1.5, 0.7), -1.5, 0.7);
}

TEST(FromCoeffs, OneDimIdentityDomain) {
    expect1DBitExact(fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0), -1.0, 1.0);
}

TEST(FromCoeffs, OneDimFusedDomain) {
    expect1DBitExact(fit<8>([](double x) { return std::cos(x); }, 0.2, 0.9, FuseAlways{}), 0.2, 0.9);
}

TEST(FromCoeffs, OneDimComplexOutput) {
    const auto fe = fit<8>([](double x) { return std::complex<double>(std::cos(x), x * x); }, -0.4, 1.3);
    const auto re = rebuild1D(fe);
    for (int i = 0; i <= 64; ++i) {
        const double x = -0.4 + 1.7 * i / 64.0;
        EXPECT_EQ(fe(x), re(x));
    }
}

TEST(FromCoeffs, NDCompileTimeCount) {
    auto f = [](const std::array<double, 2> &x) {
        return std::array<double, 2>{std::cos(x[0]) * x[1], std::exp(x[0] - x[1])};
    };
    expectNDBitExact(fit<7>(f, std::array<double, 2>{-1.0, 0.5}, std::array<double, 2>{1.0, 1.5}));
}

TEST(FromCoeffs, NDRuntimeCount) {
    auto f = [](const std::array<double, 2> &x) {
        return std::array<double, 1>{std::sin(x[0] + 2 * x[1])};
    };
    expectNDBitExact(fit(f, 6, std::array<double, 2>{-1.0, 0.5}, std::array<double, 2>{1.0, 1.5}));
}

// Positive control: the bit-exact comparison above must fire when the stored
// state differs, else the tests prove nothing.
TEST(FromCoeffs, PerturbedCoefficientIsDetected) {
    const auto fe = fit<9>([](double x) { return std::exp(x); }, 0.3, 2.1);
    using FE = decltype(fe);
    auto c = fe.coeffs();
    c[c.size() - 1] += 1e-12;
    const FE re(fromCoeffs, c, fe.domainInvSpan(), fe.domainSumEndpoints());
    EXPECT_NE(fe(1.0), re(1.0));
    const FE shifted(fromCoeffs, fe.coeffs(), fe.domainInvSpan(), fe.domainSumEndpoints() + 1e-9);
    EXPECT_NE(fe(1.0), shifted(1.0));
}

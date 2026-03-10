#include "polyfit/polyfit.hpp"
#include <cmath>
#if __cplusplus >= 202002L
#include <numbers>
#endif
#include <gtest/gtest.h>

using namespace poly_eval;

TEST(Utils, BjorckNewtonRoundtrip) {
    // Polynomial p(x) = 1 + 2x + 3x^2
    using X = double;
    using Y = double;
    const std::vector<X> nodes = {0.0, 1.0, 2.0};
    const std::vector<Y> values = {1.0, 6.0, 17.0};

    // Create Buffers (runtime-size Buffers are std::vector)
    poly_eval::Buffer<X, 0> bx_buf;
    poly_eval::Buffer<Y, 0> by_buf;
    bx_buf.assign(nodes.begin(), nodes.end());
    by_buf.assign(values.begin(), values.end());

    auto newton = detail::bjorck_pereyra<0, X, Y>(bx_buf, by_buf);
    auto mono = detail::newton_to_monomial<0, X, Y>(newton, bx_buf);

    // Evaluate monomial coefficients at original nodes and compare
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        double x = nodes[i];
        double eval = 0.0;
        double pow = 1.0;
        for (std::size_t k = 0; k < mono.size(); ++k) {
            eval += mono[k] * pow;
            pow *= x;
        }
        EXPECT_NEAR(eval, values[i], 1e-12);
    }
}

TEST(Utils, BjorckNewtonHighDegreeAccuracy) {
    // 9-coefficient polynomial p(x) = x^8 + x^7 + ... + x + 1
    // on Chebyshev nodes of the second kind in [-1, 1].
    // At 9 coefficients the monomial Vandermonde condition number is modest (~256),
    // so compensated arithmetic recovers coefficients to near machine precision.
    constexpr std::size_t nCoeffs = 9;

    std::vector<double> expectedCoeffs(nCoeffs, 1.0);

    // Generate Chebyshev nodes: x_k = cos(k * pi / (nCoeffs-1)), k = 0..nCoeffs-1
    poly_eval::Buffer<double, 0> nodes(nCoeffs);
#if __cplusplus >= 202002L
    constexpr double pi = std::numbers::pi;
#else
    constexpr double pi = 3.14159265358979323846;
#endif
    for (std::size_t k = 0; k < nCoeffs; ++k)
        nodes[k] = std::cos(static_cast<double>(k) * pi / static_cast<double>(nCoeffs - 1));

    // Evaluate via Horner: p(x) = 1 + x*(1 + x*(1 + ... + x*(1 + x)))
    poly_eval::Buffer<double, 0> values(nCoeffs);
    for (std::size_t i = 0; i < nCoeffs; ++i) {
        double x = nodes[i];
        double val = 1.0;
        for (std::size_t k = 1; k < nCoeffs; ++k) val = std::fma(val, x, 1.0);
        values[i] = val;
    }

    auto newton = detail::bjorck_pereyra<0, double, double>(nodes, values);
    auto mono = detail::newton_to_monomial<0, double, double>(newton, nodes);

    double maxErr = 0.0;
    for (std::size_t k = 0; k < nCoeffs; ++k) {
        double err = std::abs(mono[k] - expectedCoeffs[k]);
        maxErr = std::max(maxErr, err);
    }
    // Verify compensation is near machine precision
    EXPECT_LT(maxErr, 1e-13) << "max coefficient error = " << maxErr;
}

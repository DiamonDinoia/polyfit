#include <gtest/gtest.h>
#include "polyfit/fast_eval.hpp"

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

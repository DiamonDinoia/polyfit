#include <gtest/gtest.h>
#include "polyfit/fast_eval.hpp"

TEST(FuncEvalMany, PackedMatchesIndividual) {
    using namespace poly_eval;
    auto f1 = [](double x) { return std::sin(x); };
    auto f2 = [](double x) { return std::cos(2 * x); };
    double a = -1.0, b = 1.0;

    auto fe1 = make_func_eval<8>(f1, a, b);
    auto fe2 = make_func_eval<8>(f2, a, b);

    auto packed = make_func_eval_many(fe1, fe2);

    double x = 0.3;
    auto p = packed(x);
    auto r1 = fe1(x);
    auto r2 = fe2(x);
    EXPECT_NEAR(p[0], r1, 1e-6);
    EXPECT_NEAR(p[1], r2, 1e-6);

    // bulk evaluation with per-polynomial inputs
    std::array<double, 2> xs = {0.1, 0.2};
    auto p2 = packed(xs);
    EXPECT_NEAR(p2[0], fe1(xs[0]), 1e-6);
    EXPECT_NEAR(p2[1], fe2(xs[1]), 1e-6);
}

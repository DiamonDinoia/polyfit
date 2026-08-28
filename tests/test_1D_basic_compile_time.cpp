// test_1D_basic_compile_time.cpp: compile-time coefficient count fitting tests

#include "test_1D_helpers.hpp"

static std::mt19937 gen(42);

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

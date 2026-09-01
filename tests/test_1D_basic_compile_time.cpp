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

// horner on an xsimd batch returns a batch; every lane equals the scalar Horner result.
TEST(PolyEval, HornerBatchInputMatchesScalar) {
    using B = xsimd::batch<double>;
    constexpr std::array<double, 5> c{1.0, -2.0, 0.5, 3.0, -1.0};
    alignas(B) std::array<double, B::size> xs{};
    for (std::size_t i = 0; i < B::size; ++i) xs[i] = -1.0 + 2.0 * double(i) / double(B::size);
    const auto x = B::load_aligned(xs.data());
    const auto fixed = poly_eval::horner<5>(x, c.data());
    const auto dynamic = poly_eval::horner(x, c.data(), c.size());
    static_assert(std::is_same_v<decltype(fixed), const B>);
    static_assert(std::is_same_v<decltype(dynamic), const B>);
    for (std::size_t i = 0; i < B::size; ++i) {
        EXPECT_EQ(fixed.get(i), poly_eval::horner<5>(xs[i], c.data()));
        EXPECT_EQ(dynamic.get(i), poly_eval::horner(xs[i], c.data(), c.size()));
    }
}

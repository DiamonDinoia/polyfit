#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <tuple>
#include <vector>

#include "polyfit/fast_eval.hpp"

//------------------------------------------------------------------------------
// Type-parameterized test suite for FuncEvalMany with group sizes
// 2, 3, 4, 5, 12, 16. Each TEST exercises all public operators.
//------------------------------------------------------------------------------

// RNG & function helper
static std::mt19937 &rng() {
    static std::mt19937 gen(42);
    return gen;
}
static double rand_uniform() {
    static std::uniform_real_distribution<double> dist(-1.0, 1.0);
    return dist(rng());
}
static double func(double x) { return std::cos(x); }

// Helper to build a FuncEvalMany<N>
template <std::size_t N, std::size_t... Is> auto make_group_impl(std::index_sequence<Is...>) {
    const auto make_one = [] { return poly_eval::make_func_eval(func, 16, -1.0, 1.0); };
    return poly_eval::make_func_eval_many((static_cast<void>(Is), make_one())...);
}
template <std::size_t N> auto make_group() { return make_group_impl<N>(std::make_index_sequence<N>{}); }

// Test fixture: N is captured via integral_constant
template <typename IC> class FastEvalManyTest : public ::testing::Test {
  public:
    static constexpr std::size_t N = IC::value;
    using Group = decltype(make_group<N>());
    Group group = make_group<N>();
};

using GroupSizes = ::testing::Types<std::integral_constant<std::size_t, 2>, std::integral_constant<std::size_t, 3>,
                                    std::integral_constant<std::size_t, 4>, std::integral_constant<std::size_t, 5>,
                                    std::integral_constant<std::size_t, 12>, std::integral_constant<std::size_t, 16>>;
TYPED_TEST_SUITE(FastEvalManyTest, GroupSizes);

// Introspection: size() and degree()
TYPED_TEST(FastEvalManyTest, Introspection) {
    EXPECT_EQ(this->group.size(), this->N);
    EXPECT_EQ(this->group.degree(), 16u);
}

// Pointer overload, small batch (N points)
TYPED_TEST(FastEvalManyTest, PointerSmall) {
    const std::size_t M = this->N;
    std::vector<double> xs(M), out(M * M);
    for (std::size_t i = 0; i < M; ++i)
        xs[i] = rand_uniform();
    this->group(xs.data(), out.data(), M);
    stdex::mdspan<double, stdex::extents<std::size_t, stdex::dynamic_extent, M>> view(out.data(), M);
    for (std::size_t i = 0; i < M; ++i) {
        double e = func(xs[i]);
        for (std::size_t f = 0; f < M; ++f) {
            EXPECT_NEAR(view(i, f), e, 1e-15);
        }
    }
}

// Pointer overload, large batch (1000 points)
TYPED_TEST(FastEvalManyTest, PointerLarge) {
    const std::size_t M = 1000;
    const std::size_t F = this->N;
    std::vector<double> xs(M), out(M * F);
    for (std::size_t i = 0; i < M; ++i)
        xs[i] = rand_uniform();
    this->group(xs.data(), out.data(), M);
    stdex::mdspan<double, stdex::extents<std::size_t, stdex::dynamic_extent, F>> view(out.data(), M);
    for (std::size_t i = 0; i < M; ++i) {
        double e = func(xs[i]);
        for (std::size_t f = 0; f < F; ++f) {
            EXPECT_NEAR(view(i, f), e, 1e-13);
        }
    }
}

// Scalar overload
TYPED_TEST(FastEvalManyTest, Scalar) {
    double x = rand_uniform();
    auto out = this->group(x);
    for (std::size_t f = 0; f < this->N; ++f)
        EXPECT_NEAR(out[f], func(x), 1e-15);
}

// Array overload
TYPED_TEST(FastEvalManyTest, Array) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i)
        xs[i] = rand_uniform();
    auto out = this->group(xs);
    for (std::size_t f = 0; f < this->N; ++f)
        EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

// Variadic overload
TYPED_TEST(FastEvalManyTest, Variadic) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i)
        xs[i] = rand_uniform();
    auto out = std::apply([&](auto... vals) { return this->group(vals...); }, xs);
    for (std::size_t f = 0; f < this->N; ++f)
        EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

// Tuple overload
TYPED_TEST(FastEvalManyTest, Tuple) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i)
        xs[i] = rand_uniform();
    auto tup = std::apply([](auto... vals) { return std::make_tuple(vals...); }, xs);
    auto out = this->group(tup);
    for (std::size_t f = 0; f < this->N; ++f)
        EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

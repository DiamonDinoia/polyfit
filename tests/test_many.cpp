#include <array>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <random>
#include <tuple>
#include <vector>

#include "polyfit/polyfit.hpp"

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
static double mixed_func(double x) { return std::sin(x); }
static std::complex<double> complex_func(double x) { return {std::sin(x), std::cos(x)}; }

template<std::size_t I> static auto makeMixedEval() {
    if constexpr ((I % 2) == 0) {
        return poly_eval::fit(mixed_func, 16, -1.0, 1.0);
    } else {
        return poly_eval::fit(mixed_func, 16, 0.0, 1.0);
    }
}

// Helper to build a FuncEvalMany<N>
template<std::size_t N, std::size_t... Is> auto makeGroupImpl(std::index_sequence<Is...>) {
    const auto make_one = [] {
        return poly_eval::fit(func, 16, -1.0, 1.0);
    };
    return poly_eval::pack((static_cast<void>(Is), make_one())...);
}
template<std::size_t N> auto makeGroup() { return makeGroupImpl<N>(std::make_index_sequence<N>{}); }

template<std::size_t N, std::size_t... Is> auto makeMixedEvalArrayImpl(std::index_sequence<Is...>) {
    return std::array{makeMixedEval<Is>()...};
}

template<std::size_t N> auto makeMixedEvalArray() {
    return makeMixedEvalArrayImpl<N>(std::make_index_sequence<N>{});
}

template<std::size_t N> auto makeMixedGroup() {
    return std::apply([](const auto &...evals) { return poly_eval::pack(evals...); }, makeMixedEvalArray<N>());
}

template<std::size_t N, std::size_t... Is> auto makeComplexGroupImpl(std::index_sequence<Is...>) {
    const auto make_one = [] {
        return poly_eval::fit(complex_func, 16, -1.0, 1.0);
    };
    return poly_eval::pack((static_cast<void>(Is), make_one())...);
}

template<std::size_t N> auto makeComplexGroup() { return makeComplexGroupImpl<N>(std::make_index_sequence<N>{}); }

// Test fixture: N is captured via integral_constant
template<typename IC> class FastEvalManyTest : public ::testing::Test {
  public:
    static constexpr std::size_t N = IC::value;
    using Group = decltype(makeGroup<N>());
    Group evals = makeGroup<N>();
};

using GroupSizes = ::testing::Types<std::integral_constant<std::size_t, 2>, std::integral_constant<std::size_t, 3>,
                                    std::integral_constant<std::size_t, 4>, std::integral_constant<std::size_t, 5>,
                                    std::integral_constant<std::size_t, 12>, std::integral_constant<std::size_t, 16>>;
TYPED_TEST_SUITE(FastEvalManyTest, GroupSizes);

// Introspection: size() and coefficient count
TYPED_TEST(FastEvalManyTest, Introspection) {
    EXPECT_EQ(this->evals.size(), this->N);
    EXPECT_EQ(this->evals.nCoeffs(), 16u);
}

// Pointer overload, small batch (N points)
TYPED_TEST(FastEvalManyTest, PointerSmall) {
    constexpr std::size_t M = TestFixture::N;
    std::vector<double> xs(M), out(M * M);
    for (std::size_t i = 0; i < M; ++i) xs[i] = rand_uniform();
    this->evals(xs.data(), out.data(), M);
    stdex::mdspan<double, stdex::extents<std::size_t, stdex::dynamic_extent, M>> view(out.data(), M);
    for (std::size_t i = 0; i < M; ++i) {
        double e = func(xs[i]);
        for (std::size_t f = 0; f < M; ++f) {
            EXPECT_NEAR((view[std::array<std::size_t, 2>{i, f}]), e, 1e-15);
        }
    }
}

// Pointer overload, large batch (1000 points)
TYPED_TEST(FastEvalManyTest, PointerLarge) {
    const std::size_t M = 1000;
    constexpr std::size_t F = TestFixture::N;
    std::vector<double> xs(M), out(M * F);
    for (std::size_t i = 0; i < M; ++i) xs[i] = rand_uniform();
    this->evals(xs.data(), out.data(), M);
    stdex::mdspan<double, stdex::extents<std::size_t, stdex::dynamic_extent, F>> view(out.data(), M);
    for (std::size_t i = 0; i < M; ++i) {
        double e = func(xs[i]);
        for (std::size_t f = 0; f < F; ++f) {
            EXPECT_NEAR((view[std::array<std::size_t, 2>{i, f}]), e, 1e-13);
        }
    }
}

// Scalar overload
TYPED_TEST(FastEvalManyTest, Scalar) {
    double x = rand_uniform();
    auto out = this->evals(x);
    for (std::size_t f = 0; f < this->N; ++f) EXPECT_NEAR(out[f], func(x), 1e-15);
}

// Array overload
TYPED_TEST(FastEvalManyTest, Array) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i) xs[i] = rand_uniform();
    auto out = this->evals(xs);
    for (std::size_t f = 0; f < this->N; ++f) EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

// Variadic overload
TYPED_TEST(FastEvalManyTest, Variadic) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i) xs[i] = rand_uniform();
    auto out = std::apply([&](auto... vals) { return this->evals(vals...); }, xs);
    for (std::size_t f = 0; f < this->N; ++f) EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

// Tuple overload
TYPED_TEST(FastEvalManyTest, Tuple) {
    std::array<double, TestFixture::N> xs;
    for (std::size_t i = 0; i < TestFixture::N; ++i) xs[i] = rand_uniform();
    auto tup = std::apply([](auto... vals) { return std::make_tuple(vals...); }, xs);
    auto out = this->evals(tup);
    for (std::size_t f = 0; f < this->N; ++f) EXPECT_NEAR(out[f], func(xs[f]), 1e-15);
}

// ----- Public API shape assertions -----

TEST(FuncEvalManyAPI, ExposesOnlyUsefulCount) {
    using Group2 = decltype(makeGroup<2>());
    using Group5 = decltype(makeGroup<5>());

    EXPECT_EQ(Group2::COUNT, 2u);
    EXPECT_EQ(Group5::COUNT, 5u);
}

// ----- FuncEvalMany truncation tests -----

TEST(FuncEvalManyTruncate, BasicTruncation) {
    auto group = makeGroup<4>();
    EXPECT_EQ(group.nCoeffs(), 16u);

    // Truncate with a generous threshold — should reduce the coefficient count
    group.truncate(1e-8);
    EXPECT_LE(group.nCoeffs(), 16u);

    // Verify accuracy is maintained
    double x = 0.3;
    auto out = group(x);
    double expected = func(x);
    for (std::size_t f = 0; f < 4; ++f) EXPECT_NEAR(out[f], expected, 1e-7);
}

TEST(FuncEvalManyTruncate, PreservesAccuracy) {
    auto group = makeGroup<3>();

    // Truncate with a tight threshold
    group.truncate(1e-14);

    // Verify across multiple points
    for (int i = 0; i < 50; ++i) {
        double x = rand_uniform();
        auto out = group(x);
        for (std::size_t f = 0; f < 3; ++f) EXPECT_NEAR(out[f], func(x), 1e-13);
    }
}

// Non-SIMD-multiple sizes should still work correctly with padding
TEST(FuncEvalManySIMD, NonMultipleSizes) {
    // Size 3 and 5 are not multiples of typical SIMD widths (2 or 4)
    auto group3 = makeGroup<3>();
    auto group5 = makeGroup<5>();

    double x = 0.42;
    auto out3 = group3(x);
    auto out5 = group5(x);
    double expected = func(x);

    for (std::size_t f = 0; f < 3; ++f) EXPECT_NEAR(out3[f], expected, 1e-14);
    for (std::size_t f = 0; f < 5; ++f) EXPECT_NEAR(out5[f], expected, 1e-14);
}

TEST(FuncEvalManyMixedDomain, SameInputMatchesIndividual) {
    constexpr std::size_t count = 4;
    auto group = makeMixedGroup<count>();
    auto evals = makeMixedEvalArray<count>();

    const double x = 0.25;
    const auto out = group(x);
    for (std::size_t i = 0; i < count; ++i) {
        EXPECT_NEAR(out[i], evals[i](x), 1e-15);
    }
}

TEST(FuncEvalManyMixedDomain, BatchMatchesIndividual) {
    constexpr std::size_t count = 4;
    const std::size_t num_points = 32;
    auto group = makeMixedGroup<count>();
    auto evals = makeMixedEvalArray<count>();

    std::vector<double> xs(num_points), out(num_points * count);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (double &x : xs) x = dist(rng());

    group(xs.data(), out.data(), num_points);
    stdex::mdspan<double, stdex::extents<std::size_t, stdex::dynamic_extent, count>> view(out.data(), num_points);
    for (std::size_t p = 0; p < num_points; ++p) {
        for (std::size_t i = 0; i < count; ++i) {
            EXPECT_NEAR((view[std::array<std::size_t, 2>{p, i}]), evals[i](xs[p]), 1e-13);
        }
    }
}

TEST(FuncEvalManyTypeMismatch, ComplexOutputBulkMatchesScalar) {
    constexpr std::size_t count = 4;
    const std::size_t num_points = 32;
    auto group = makeComplexGroup<count>();

    std::vector<double> xs(num_points);
    std::vector<std::complex<double>> out(num_points * count);
    for (double &x : xs) x = rand_uniform();

    group(xs.data(), out.data(), num_points);
    stdex::mdspan<std::complex<double>, stdex::extents<std::size_t, stdex::dynamic_extent, count>> view(out.data(), num_points);
    for (std::size_t p = 0; p < num_points; ++p) {
        const auto expected = complex_func(xs[p]);
        for (std::size_t i = 0; i < count; ++i) {
            EXPECT_NEAR((view[std::array<std::size_t, 2>{p, i}].real()), expected.real(), 1e-13);
            EXPECT_NEAR((view[std::array<std::size_t, 2>{p, i}].imag()), expected.imag(), 1e-13);
        }
    }
}

TEST(FuncEvalManyTruncate, ComplexOutputCompilesAndPreservesAccuracy) {
    auto group = makeComplexGroup<4>();
    group.truncate(1e-12);

    const double x = 0.3;
    const auto out = group(x);
    const auto expected = complex_func(x);
    for (const auto &value : out) {
        EXPECT_NEAR(value.real(), expected.real(), 1e-12);
        EXPECT_NEAR(value.imag(), expected.imag(), 1e-12);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

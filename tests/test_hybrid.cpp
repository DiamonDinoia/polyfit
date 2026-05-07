// test_hybrid.cpp
// -----------------------------------------------------------------------------
// Hybrid Estrin/Horner evaluator tests, split out of test_horner.cpp so MSVC
// stays under its template-instantiation-context limit. The compile-time
// `hybrid<D>` and `hybrid_transposed<N>` paths each contain several layers of
// nested `poet::static_for`, and instantiating them for D = 1..32 alongside
// the rest of test_horner.cpp tripped MSVC's C1202 limit.
// -----------------------------------------------------------------------------

#include <array>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

#include "polyfit/polyfit.hpp"

#include <poet/poet.hpp>

using std::size_t;

namespace {

template<typename T> constexpr T eps = std::numeric_limits<T>::epsilon() * T(100);

std::mt19937 rng(42);
std::uniform_real_distribution<double> uni_dist(-1.0, 1.0);

template<typename T> T uni() { return static_cast<T>(uni_dist(rng)); }

template<typename T> std::vector<T> random_vector(size_t n) {
    std::vector<T> v(n);
    for (auto &x : v) x = uni<T>();
    return v;
}

template<typename Output, typename Input>
constexpr Output naive_horner_scalar(Input x, const Output *c, size_t n) {
    Output acc = c[0];
    for (size_t i = 1; i < n; ++i) acc = acc * x + c[i];
    return acc;
}

template<typename T> class HybridTyped : public testing::Test {};
typedef testing::Types<float, double> FloatingTypes;
TYPED_TEST_SUITE(HybridTyped, FloatingTypes);

} // namespace

// Hybrid (mixed Estrin/Horner) — runtime degree dispatches in [1, 32], falls
// back to serial Horner above 32. Verify numerical equivalence with the
// compile-time path across the dispatch boundary.
TYPED_TEST(HybridTyped, RuntimeDispatch_MatchesCompileTime) {
    using T = TypeParam;
    for (size_t N : {size_t(1), size_t(7), size_t(16), size_t(32), size_t(33)}) {
        std::vector<T> c = random_vector<T>(N);
        T x = uni<T>();
        T rt = poly_eval::hybrid(x, c.data(), c.size());
        T ex = naive_horner_scalar(x, c.data(), c.size());
        EXPECT_NEAR(rt, ex, eps<T> * static_cast<T>(N));
    }
}

TYPED_TEST(HybridTyped, CompileTime_MatchesNaive) {
    using T = TypeParam;
    auto check = [&](auto D) {
        std::vector<T> c = random_vector<T>(D);
        T x = uni<T>();
        T ct = poly_eval::hybrid<D>(x, c.data());
        T ex = naive_horner_scalar(x, c.data(), D);
        EXPECT_NEAR(ct, ex, eps<T> * static_cast<T>(D));
    };
    poet::static_for<1, 17>(check);
    poet::static_for<17, 33>(check);
}

// Hybrid transposed (SIMD over a precomputed transposed coefficient buffer)
TYPED_TEST(HybridTyped, Transposed_MatchesNaive) {
    using T = TypeParam;
    using Batch = xsimd::batch<T>;
    auto check = [&](auto D) {
        constexpr size_t N = D;
        std::vector<T> c = random_vector<T>(N);
        constexpr size_t Sz = poly_eval::hybrid_transposed_size<N, Batch>();
        alignas(Batch::arch_type::alignment()) T c_trans[Sz];
        poly_eval::hybrid_transpose_coeffs<N, Batch>(c.data(), c_trans);
        T x = uni<T>();
        T got = poly_eval::hybrid_transposed<N, Batch>(x, c_trans);
        T ex = naive_horner_scalar(x, c.data(), N);
        EXPECT_NEAR(got, ex, eps<T> * static_cast<T>(N));
    };
    poet::static_for<5, 19>(check);
    poet::static_for<19, 33>(check);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

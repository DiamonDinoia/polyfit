// test_block_size.cpp
// -----------------------------------------------------------------------------
// Tests for the parametric Hybrid block-size knob:
//   - `HybridK<K>` tag override on `hybrid_impl_ct` / `hybrid` /
//     `hybrid_transposed_size` / `hybrid_transposed` produces numerical
//     parity with the default (heuristic) path.
//   - `optimal_block_size<NCOEFFS, SIMD_W, NREG, EvalPolicy>()` returns
//     values in the expected range and respects the documented policy
//     ordering (Throughput K >= Latency K when SIMD_W >= 2, etc.).
// -----------------------------------------------------------------------------

#include <array>
#include <cstddef>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <vector>

#include "polyfit/polyfit.hpp"

#include <poet/poet.hpp>

using std::size_t;

namespace {

template<typename T> constexpr T eps = std::numeric_limits<T>::epsilon() * T(100);

std::mt19937 rng(1729);
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

} // namespace

// -----------------------------------------------------------------------------
// optimal_block_size: compile-time shape sanity for representative microarchs.
// -----------------------------------------------------------------------------
namespace {
using poly_eval::EvalPolicy;
using poly_eval::detail::optimal_block_size;

// AVX-512: SIMD_W = 8, NREG = 32
constexpr size_t kAvx512W = 8;
constexpr size_t kAvx512R = 32;
static_assert(optimal_block_size<3, kAvx512W, kAvx512R, EvalPolicy::Latency>() == 3, "");
static_assert(optimal_block_size<4, kAvx512W, kAvx512R, EvalPolicy::Balanced>() == 4, "");
static_assert(optimal_block_size<5, kAvx512W, kAvx512R, EvalPolicy::Throughput>() == 4, "");
static_assert(optimal_block_size<16, kAvx512W, kAvx512R, EvalPolicy::Throughput>() == 15, "");
static_assert(optimal_block_size<16, kAvx512W, kAvx512R, EvalPolicy::Latency>() >= 2, "");
static_assert(optimal_block_size<16, kAvx512W, kAvx512R, EvalPolicy::Latency>() <= 15, "");
static_assert(optimal_block_size<16, kAvx512W, kAvx512R, EvalPolicy::Throughput>()
                  >= optimal_block_size<16, kAvx512W, kAvx512R, EvalPolicy::Latency>(),
              "");

// AVX2: SIMD_W = 4, NREG = 16
constexpr size_t kAvx2W = 4;
constexpr size_t kAvx2R = 16;
static_assert(optimal_block_size<16, kAvx2W, kAvx2R, EvalPolicy::Balanced>() >= 2, "");
static_assert(optimal_block_size<16, kAvx2W, kAvx2R, EvalPolicy::Balanced>() <= 15, "");
static_assert(optimal_block_size<16, kAvx2W, kAvx2R, EvalPolicy::Throughput>() == 15, "");

// NEON-fp64: SIMD_W = 2, NREG = 32
static_assert(optimal_block_size<16, 2, 32, EvalPolicy::Throughput>() == 15, "");

// Scalar: SIMD_W = 1, any NREG; Throughput falls back to Latency.
static_assert(optimal_block_size<16, 1, 16, EvalPolicy::Throughput>()
                  == optimal_block_size<16, 1, 16, EvalPolicy::Latency>(),
              "");
static_assert(optimal_block_size<32, 1, 32, EvalPolicy::Throughput>()
                  == optimal_block_size<32, 1, 32, EvalPolicy::Latency>(),
              "");

} // namespace

TEST(OptimalBlockSize, RuntimeSanityWrapsConstexpr) {
    // The static_asserts above guard these values at compile time. Repeat the
    // key ones as runtime checks.
    constexpr size_t k8 = optimal_block_size<8, 8, 32, EvalPolicy::Throughput>();
    constexpr size_t k16 = optimal_block_size<16, 8, 32, EvalPolicy::Throughput>();
    EXPECT_EQ(k8, size_t(7));
    EXPECT_EQ(k16, size_t(15));
    const size_t latencyK = optimal_block_size<16, 8, 32, EvalPolicy::Latency>();
    EXPECT_GE(latencyK, size_t(2));
    EXPECT_LE(latencyK, size_t(15));
}

TEST(EvalPolicy, EnumRoundtrip) {
    EvalPolicy p = EvalPolicy::Throughput;
    EXPECT_EQ(static_cast<int>(p), static_cast<int>(EvalPolicy::Throughput));
    EXPECT_NE(static_cast<int>(EvalPolicy::Latency), static_cast<int>(EvalPolicy::Throughput));
    EXPECT_NE(static_cast<int>(EvalPolicy::Balanced), static_cast<int>(EvalPolicy::Latency));
}

// -----------------------------------------------------------------------------
// HybridK override: numerical parity with the default heuristic.
// -----------------------------------------------------------------------------
template<typename T> class HybridKTyped : public testing::Test {};
using FloatingTypes = testing::Types<float, double>;
TYPED_TEST_SUITE(HybridKTyped, FloatingTypes);

namespace {
template<size_t N, size_t K, typename T>
void check_hybrid_K_matches_default(const std::vector<T> &c, T x) {
    // FMA grouping changes with K, so bit-exact equality is not promised. The
    // `eps * N` tolerance of the rest of the hybrid suite is the binding criterion.
    T withK = poly_eval::detail::hybrid_impl_ct<N, T, T, T, K>(x, c.data());
    T naive = naive_horner_scalar(x, c.data(), N);
    EXPECT_NEAR(withK, naive, eps<T> * static_cast<T>(N))
        << "N=" << N << " K=" << K;
}
} // namespace

TYPED_TEST(HybridKTyped, ScalarOverride_MatchesNaive) {
    using T = TypeParam;
    auto sweep = [&](auto Ntag) {
        constexpr size_t N = Ntag;
        std::vector<T> c = random_vector<T>(N);
        T x = uni<T>();
        // K in {2, 3, 4, N-1}
        check_hybrid_K_matches_default<N, 2>(c, x);
        check_hybrid_K_matches_default<N, 3>(c, x);
        check_hybrid_K_matches_default<N, 4>(c, x);
        check_hybrid_K_matches_default<N, (N >= 3 ? N - 1 : 2)>(c, x);
    };
    sweep(std::integral_constant<size_t, 5>{});
    sweep(std::integral_constant<size_t, 6>{});
    sweep(std::integral_constant<size_t, 8>{});
    sweep(std::integral_constant<size_t, 10>{});
    sweep(std::integral_constant<size_t, 16>{});
}

TYPED_TEST(HybridKTyped, TransposedOverride_MatchesNaive) {
    using T = TypeParam;
    using Batch = xsimd::batch<T>;
    auto sweep = [&](auto Ntag) {
        constexpr size_t N = Ntag;
        constexpr size_t K = 3; // pick a non-default K; layout invariants
                                // make it valid for every N > 4.
        std::vector<T> c = random_vector<T>(N);
        constexpr size_t Sz = poly_eval::hybrid_transposed_size<N, Batch, K>();
        std::vector<T> c_trans_storage(Sz + Batch::size, T(0));
        // Align by hand; a vector is heap-aligned but not necessarily SIMD-aligned.
        const auto raw = reinterpret_cast<std::uintptr_t>(c_trans_storage.data());
        const auto algn = Batch::arch_type::alignment();
        const auto off = (algn - (raw % algn)) % algn;
        T *c_trans = reinterpret_cast<T *>(raw + off);
        poly_eval::hybrid_transpose_coeffs<N, Batch, K>(c.data(), c_trans);
        T x = uni<T>();
        T got = poly_eval::hybrid_transposed<N, Batch, T, K>(x, c_trans);
        T ex = naive_horner_scalar(x, c.data(), N);
        EXPECT_NEAR(got, ex, eps<T> * static_cast<T>(N)) << "N=" << N << " K=" << K;
    };
    sweep(std::integral_constant<size_t, 5>{});
    sweep(std::integral_constant<size_t, 8>{});
    sweep(std::integral_constant<size_t, 10>{});
    sweep(std::integral_constant<size_t, 16>{});
}

// -----------------------------------------------------------------------------
// HybridK tag plumbed through `fit(...)` and FuncEval: end-to-end parity.
// -----------------------------------------------------------------------------
TEST(HybridKTag, FitForwardsHybridK) {
    auto F = [](double x) { return std::cos(2.0 * x) + 0.3 * x * x * x; };
    constexpr int N = 12;
    auto def = poly_eval::fit<N>(F, -1.0, 1.0);
    auto pinned = poly_eval::fit<N>(F, -1.0, 1.0, poly_eval::HybridK<3>{});
    // Coefficients are basis-only, so both fits compute the same monomial coefficients.
    for (size_t i = 0; i < static_cast<size_t>(N); ++i) {
        EXPECT_EQ(def.coeffs()[i], pinned.coeffs()[i]);
    }
    // Evaluations may differ in FMA grouping but must agree to fp tolerance.
    for (double x = -1.0; x <= 1.0; x += 0.1) {
        EXPECT_NEAR(def(x), pinned(x), 1e-12);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

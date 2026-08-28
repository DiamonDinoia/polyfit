// test_hybrid_transposed.cpp: SIMD transposed-coefficient hybrid vs naive.

#include "test_hybrid_shared.hpp"

using namespace test_hybrid_detail;

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

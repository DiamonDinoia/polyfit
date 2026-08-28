// test_hybrid_runtime_dispatch.cpp: runtime degree dispatch vs compile-time.

#include "test_hybrid_shared.hpp"

using namespace test_hybrid_detail;

// Hybrid (mixed Estrin/Horner): a runtime coefficient count dispatches in
// [1, 32] and falls back to serial Horner above. Verify numerical equivalence
// with the compile-time path across the dispatch boundary.
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

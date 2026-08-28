// test_hybrid_compile_time.cpp: compile-time hybrid degrees 1..32 vs naive.

#include "test_hybrid_shared.hpp"

using namespace test_hybrid_detail;

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

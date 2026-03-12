#include <gtest/gtest.h>

#include "polyfit/internal/macros.h"
#include "polyfit/internal/numeric_utils.h"

namespace {

PF_CXX20_CONSTEVAL int ct() { return 7; }
PF_CXX20_CONSTEXPR int cx() { return 5; }

constexpr int ctv = ct();

constexpr int branch() {
    PF_IF_CONSTEVAL { return 1; }
    return 2;
}

constexpr int local() {
    PF_STATIC_CONSTEXPR_LOCAL int v = 11;
    return v;
}

static_assert(ctv == 7, "consteval/constexpr compatibility");
#if PF_HAS_CXX20
static_assert(cx() == 5, "constexpr compatibility");
#endif
static_assert(local() == 11, "local static constexpr compatibility");
static_assert(!PF_HAS_CXX23 || PF_HAS_CXX20, "C++23 implies C++20");
static_assert(poly_eval::detail::math::fma(2.0, 3.0, 1.0) == 7.0, "math::fma constexpr");
static_assert(poly_eval::detail::math::sqrt(4.0) == 2.0, "math::sqrt constexpr");

#if PF_HAS_CXX20 || (defined(_MSC_VER) && _MSC_VER >= 1925) || PF_HAS_BUILTIN(__builtin_is_constant_evaluated)
static_assert(branch() == 1, "constant-evaluated branch");
#endif

TEST(FeatureMacros, Runtime) {
    EXPECT_EQ(branch(), 2);
    EXPECT_EQ(local(), 11);
    EXPECT_EQ(ctv, 7);
    EXPECT_EQ(cx(), 5);
}

} // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

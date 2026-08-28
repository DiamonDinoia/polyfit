// Gate: under sanitizers this 16^3 = 4096-term instantiation (about 15 MiB per
// tensor coefficient) exceeds the 16 GiB CI runners. Lighter monomial TUs keep
// ND coverage under sanitizers; full degrees run in the non-sanitizer jobs.
#ifdef POLYFIT_TESTS_REDUCE_ND
#include <gtest/gtest.h>
TEST(Eval, In3Out2Deg16) { GTEST_SKIP() << "tensor size 4096 exceeds the sanitizer build budget"; }
#else
#include "test_ND_monomial.hpp"
TEST(Eval, In3Out2Deg16) { runMonomialTest<3, 2, 16>(); }
#endif

// Gate: under sanitizers this 8^3 = 512-term instantiation (about 15 MiB per
// tensor coefficient, superlinear in wall time) exceeds the 16 GiB CI runners.
// Lighter monomial TUs keep ND coverage under sanitizers; full degrees run in
// the non-sanitizer jobs.
#ifdef POLYFIT_TESTS_REDUCE_ND
#include <gtest/gtest.h>
TEST(Eval, In3Out4Deg8) { GTEST_SKIP() << "tensor size 512 exceeds the sanitizer build budget"; }
#else
#include "test_ND_monomial.hpp"
TEST(Eval, In3Out4Deg8) { runMonomialTest<3, 4, 8>(); }
#endif

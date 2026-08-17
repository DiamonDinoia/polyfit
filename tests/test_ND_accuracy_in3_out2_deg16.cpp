// Gate: this eval instantiation unrolls 16^3 = 4096 tensor terms; with sanitizers
// it needs far more memory and time than a 16 GiB CI runner allows (~15 MiB
// per tensor coefficient under ASan+UBSan, superlinear in wall time). The
// lighter mono TUs still cover the ND fit/eval machinery under sanitizers;
// full-degree coverage remains in the non-sanitizer jobs.
#ifdef POLYFIT_TESTS_REDUCE_ND
#include <gtest/gtest.h>
TEST(Eval, In3Out2Deg16) { GTEST_SKIP() << "tensor size 4096 exceeds the sanitizer build budget"; }
#else
#include "test_ND_monomial.hpp"
TEST(Eval, In3Out2Deg16) { runMonomialTest<3, 2, 16>(); }
#endif

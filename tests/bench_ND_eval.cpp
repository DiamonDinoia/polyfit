// bench_ND_eval.cpp: single-point evaluation benchmarks for std::array points.

#include "bench_ND_shared.hpp"

void runBenchNdEval(ankerl::nanobench::Bench &bench) {
    benchEvalSingle<2, 2, 16>("F:ℝ²→ℝ², D=16", bench);
    benchEvalSingle<2, 3, 16>("F:ℝ²→ℝ³, D=16", bench);
    benchEvalSingle<3, 3, 8>("F:ℝ³→ℝ³, D=8", bench);
    benchEvalSingle<3, 4, 8>("F:ℝ³→ℝ⁴, D=8", bench);
#ifndef POLYFIT_TESTS_REDUCE_ND
    // 4096-coefficient instantiations: hours under sanitizers and tens of
    // minutes under clang/MSVC optimizers; gcc Release benches keep them.
    benchEvalSingle<3, 2, 16>("F:ℝ³→ℝ², D=16", bench);
    benchEvalSingle<4, 3, 8>("F:ℝ⁴→ℝ³, D=8", bench);
    benchEvalSingle<4, 4, 8>("F:ℝ⁴→ℝ⁴, D=8", bench);
#endif
}

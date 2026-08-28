// bench_ND_misc.cpp: generic point, variadic and batch evaluation benchmarks.

#include "bench_ND_shared.hpp"

void runBenchNdMisc(ankerl::nanobench::Bench &bench) {
    benchEvalGenericPoint<2, 2, 16>("F:ℝ²→ℝ², D=16", bench);
    benchEvalVariadic2D<2, 16>("F:ℝ²→ℝ², D=16", bench);
    benchEvalBatch<2, 2, 16>("F:ℝ²→ℝ², D=16", bench);
    benchEvalGenericBatch<2, 2, 16>("F:ℝ²→ℝ², D=16", bench);
}

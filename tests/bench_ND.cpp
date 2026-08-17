// bench_ND.cpp — runner main(); benchmark groups live in bench_ND_build.cpp,
// bench_ND_eval.cpp and bench_ND_misc.cpp (one heavyweight instantiation set
// per translation unit).

#include "bench_ND_shared.hpp"

int main() {
    ankerl::nanobench::Bench benchBuildObj;
    benchBuildObj.title("fitting time").unit("build").minEpochTime(std::chrono::milliseconds(10)).batch(1);
    runBenchNdBuild(benchBuildObj);

    ankerl::nanobench::Bench benchEvalObj;
    benchEvalObj.title("throughput (single x)")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(20))
        .batch(1);
    runBenchNdEval(benchEvalObj);
    runBenchNdMisc(benchEvalObj);

    return 0;
}

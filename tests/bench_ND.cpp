#include <array>
#include <cmath>
#include <iostream>
#include <random>

#include <nanobench.h>
#include "polyfit/fast_eval.hpp"

// -----------------------------------------------------------------------------
// Mathematical setup:
//   F : ℝ^InDim → ℝ^OutDim is our analytic function (sum of cosines).
//   We build a degree-D interpolant \tilde F on [-1,1]^InDim.
//   We then pick a single random x ∈ [-1,1]^InDim; nanobench will repeatedly
//   call \tilde F(x) to measure throughput.
// -----------------------------------------------------------------------------

constexpr std::mt19937::result_type kSeed = 42;
static std::mt19937               gen{kSeed};
static std::uniform_real_distribution<double> dist(-1.0, 1.0);

// sum of cosines analytic function (same as in your tests)
template <typename Array, typename OutPut = Array>
constexpr auto sumCos(const Array &x) {
    double s = 0.0;
    for (double xi : x) s += std::cos(xi);
    OutPut y{};
    y.fill(s);
    for (std::size_t i = 1; i < y.size(); ++i)
        y[i] += y[i - 1];
    return y;
}

// -----------------------------------------------------------------------------
// Benchmark helper: throughput of evaluating the interpolant at one point
// -----------------------------------------------------------------------------
template <std::size_t InDim, std::size_t OutDim, std::size_t Degree>
void benchEvalSingle(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input  = std::array<double, InDim>;
    using Output = std::array<double, OutDim>;

    // 1) pick one random x in [-1,1]^InDim
    Input x{};
    for (auto &xi : x) xi = dist(gen);

    // 2) build interpolant once
    Input a{}; a.fill(-1.0);
    Input b{}; b.fill( 1.0);
    auto approx = poly_eval::make_func_eval<Degree>(sumCos<Input, Output>, a, b);

    // 3) let nanobench call approx(x) repeatedly
    bench.run(label, [&] {
        ankerl::nanobench::doNotOptimizeAway(approx(x));
    });
}

// -----------------------------------------------------------------------------
// Benchmark helper: time to construct (fit) the interpolant
// -----------------------------------------------------------------------------
template <std::size_t InDim, std::size_t OutDim, std::size_t Degree>
void benchBuildSingle(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input  = std::array<double, InDim>;
    using Output = std::array<double, OutDim>;

    Input a{}; a.fill(-1.0);
    Input b{}; b.fill( 1.0);

    bench.run(label + " build", [&] {
        ankerl::nanobench::doNotOptimizeAway(poly_eval::make_func_eval<Degree>(sumCos<Input, Output>, a, b));
    });
}

// -----------------------------------------------------------------------------
// main: configure nanobench and run all cases on a single x
// -----------------------------------------------------------------------------
int main() {
    ankerl::nanobench::Bench benchEvalObj;

    // — Fitting (build) time on the same single-run basis —
    ankerl::nanobench::Bench benchBuildObj;
    benchBuildObj.title("fitting time")
                 .unit("build")
                 .minEpochTime(std::chrono::milliseconds(10))
                 .batch(1);

    benchBuildSingle<2,2,16>("F:ℝ²→ℝ², D=16", benchBuildObj);
    benchBuildSingle<2,3,16>("F:ℝ²→ℝ³, D=16", benchBuildObj);
    benchBuildSingle<3,2,16>("F:ℝ³→ℝ², D=16", benchBuildObj);
    benchBuildSingle<3,3, 8>("F:ℝ³→ℝ³, D=8" , benchBuildObj);
    benchBuildSingle<3,4, 8>("F:ℝ³→ℝ⁴, D=8" , benchBuildObj);
    benchBuildSingle<4,3, 8>("F:ℝ⁴→ℝ³, D=8" , benchBuildObj);
    benchBuildSingle<4,4, 8>("F:ℝ⁴→ℝ⁴, D=8" , benchBuildObj);

    // — Evaluation throughput on one sample x —
    benchEvalObj.title("throughput (single x)")
                 .unit("eval")
                 .minEpochTime(std::chrono::milliseconds(20))
                 .batch(1);  // each run is a single evaluation

    benchEvalSingle<2,2,16>("F:ℝ²→ℝ², D=16", benchEvalObj);
    benchEvalSingle<2,3,16>("F:ℝ²→ℝ³, D=16", benchEvalObj);
    benchEvalSingle<3,2,16>("F:ℝ³→ℝ², D=16", benchEvalObj);
    benchEvalSingle<3,3, 8>("F:ℝ³→ℝ³, D=8" , benchEvalObj);
    benchEvalSingle<3,4, 8>("F:ℝ³→ℝ⁴, D=8" , benchEvalObj);
    benchEvalSingle<4,3, 8>("F:ℝ⁴→ℝ³, D=8" , benchEvalObj);
    benchEvalSingle<4,4, 8>("F:ℝ⁴→ℝ⁴, D=8" , benchEvalObj);


    return 0;
}

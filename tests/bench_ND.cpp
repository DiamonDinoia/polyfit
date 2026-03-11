#include <array>
#include <cmath>
#include <iostream>
#include <random>

#include "polyfit/polyfit.hpp"
#include <nanobench.h>

constexpr std::mt19937::result_type kSeed = 42;
static std::mt19937 gen{kSeed};
static std::uniform_real_distribution<double> dist(-1.0, 1.0);

template<typename Array, typename OutPut = Array> constexpr auto sumCos(const Array &x) {
    double s = 0.0;
    for (double xi : x) s += std::cos(xi);
    OutPut y{};
    y.fill(s);
    for (std::size_t i = 1; i < y.size(); ++i) y[i] += y[i - 1];
    return y;
}

template<class Input> constexpr Input unitCubeMin() {
    Input a{};
    a.fill(-1.0);
    return a;
}

template<class Input> constexpr Input unitCubeMax() {
    Input b{};
    b.fill(1.0);
    return b;
}

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalSingle(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;

    Input x{};
    for (auto &xi : x) xi = dist(gen);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label, [&] { ankerl::nanobench::doNotOptimizeAway(approx(x)); });
}

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchBuildSingle(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;

    bench.run(label + " build", [&] {
        ankerl::nanobench::doNotOptimizeAway(
            poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>()));
    });
}

int main() {
    ankerl::nanobench::Bench benchEvalObj;
    ankerl::nanobench::Bench benchBuildObj;
    benchBuildObj.title("fitting time").unit("build").minEpochTime(std::chrono::milliseconds(10)).batch(1);

    benchBuildSingle<2, 2, 16>("F:ℝ²→ℝ², D=16", benchBuildObj);
    benchBuildSingle<2, 3, 16>("F:ℝ²→ℝ³, D=16", benchBuildObj);
    benchBuildSingle<3, 2, 16>("F:ℝ³→ℝ², D=16", benchBuildObj);
    benchBuildSingle<3, 3, 8>("F:ℝ³→ℝ³, D=8", benchBuildObj);
    benchBuildSingle<3, 4, 8>("F:ℝ³→ℝ⁴, D=8", benchBuildObj);
    benchBuildSingle<4, 3, 8>("F:ℝ⁴→ℝ³, D=8", benchBuildObj);
    benchBuildSingle<4, 4, 8>("F:ℝ⁴→ℝ⁴, D=8", benchBuildObj);

    benchEvalObj.title("throughput (single x)")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(20))
        .batch(1);

    benchEvalSingle<2, 2, 16>("F:ℝ²→ℝ², D=16", benchEvalObj);
    benchEvalSingle<2, 3, 16>("F:ℝ²→ℝ³, D=16", benchEvalObj);
    benchEvalSingle<3, 2, 16>("F:ℝ³→ℝ², D=16", benchEvalObj);
    benchEvalSingle<3, 3, 8>("F:ℝ³→ℝ³, D=8", benchEvalObj);
    benchEvalSingle<3, 4, 8>("F:ℝ³→ℝ⁴, D=8", benchEvalObj);
    benchEvalSingle<4, 3, 8>("F:ℝ⁴→ℝ³, D=8", benchEvalObj);
    benchEvalSingle<4, 4, 8>("F:ℝ⁴→ℝ⁴, D=8", benchEvalObj);

    return 0;
}

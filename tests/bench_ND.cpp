#include <cmath>
#include <iostream>
#include <random>

#include "test_ND_shared.hpp"
#include <nanobench.h>

constexpr std::mt19937::result_type kSeed = 42;
static std::mt19937 gen{kSeed};
static std::uniform_real_distribution<double> dist(-1.0, 1.0);

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
void benchEvalGenericPoint(const std::string &label, ankerl::nanobench::Bench &bench) {
    using CanonicalInput = std::array<double, IN_DIM>;
    using Input = FixedVec<double, IN_DIM>;
    using Output = FixedVec<double, OUT_DIM>;

    Input x{};
    for (std::size_t i = 0; i < IN_DIM; ++i) x[i] = dist(gen);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label + " generic point", [&] { ankerl::nanobench::doNotOptimizeAway(approx(x)); });

    CanonicalInput canonical{};
    for (std::size_t i = 0; i < IN_DIM; ++i) canonical[i] = x[i];
    bench.run(label + " canonical point", [&] { ankerl::nanobench::doNotOptimizeAway(approx(canonical)); });
}

template<std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalVariadic2D(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = std::array<double, 2>;
    using Output = std::array<double, OUT_DIM>;

    Input x{};
    for (auto &xi : x) xi = dist(gen);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label + " array", [&] { ankerl::nanobench::doNotOptimizeAway(approx(x)); });
    bench.run(label + " variadic", [&] { ankerl::nanobench::doNotOptimizeAway(approx(x[0], x[1])); });
}

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalBatch(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;

    constexpr std::size_t count = 1024;
    std::vector<Input> pts(count);
    for (auto &pt : pts)
        for (auto &xi : pt) xi = dist(gen);
    std::vector<Output> out(count);
    std::vector<Output> baseline(count);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label + " manual loop", [&] {
        for (std::size_t i = 0; i < pts.size(); ++i) baseline[i] = approx(pts[i]);
        ankerl::nanobench::doNotOptimizeAway(baseline.data());
    });

    bench.run(label + " canonical batch", [&] {
        approx(pts.data(), out.data(), pts.size());
        ankerl::nanobench::doNotOptimizeAway(out.data());
    });

#if defined(__cpp_lib_span) && (__cpp_lib_span >= 202002L)
    bench.run(label + " span batch", [&] {
        approx(std::span<const Input>(pts), std::span<Output>(out));
        ankerl::nanobench::doNotOptimizeAway(out.data());
    });
#endif
}

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalGenericBatch(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = FixedVec<double, IN_DIM>;
    using Output = FixedVec<double, OUT_DIM>;

    constexpr std::size_t count = 1024;
    std::vector<Input> pts(count);
    for (auto &pt : pts)
        for (std::size_t i = 0; i < IN_DIM; ++i) pt[i] = dist(gen);
    std::vector<Output> out(count);
    std::vector<Output> baseline(count);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label + " generic loop", [&] {
        for (std::size_t i = 0; i < pts.size(); ++i) baseline[i] = approx(pts[i]);
        ankerl::nanobench::doNotOptimizeAway(baseline.data());
    });

    bench.run(label + " container batch", [&] {
        approx(pts, out);
        ankerl::nanobench::doNotOptimizeAway(out.data());
    });
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
    benchEvalGenericPoint<2, 2, 16>("F:ℝ²→ℝ², D=16", benchEvalObj);
    benchEvalVariadic2D<2, 16>("F:ℝ²→ℝ², D=16", benchEvalObj);
    benchEvalBatch<2, 2, 16>("F:ℝ²→ℝ², D=16", benchEvalObj);
    benchEvalGenericBatch<2, 2, 16>("F:ℝ²→ℝ², D=16", benchEvalObj);

    return 0;
}

#pragma once

// Shared helpers for the bench_ND translation units. Each ND fit is a
// heavyweight template instantiation, so the benchmark groups live one per
// .cpp; the runner main() lives in bench_ND.cpp.

#include <array>
#include <cmath>
#include <chrono>
#include <random>
#include <string>
#include <vector>

#include "test_ND_shared.hpp"
#include <nanobench.h>

inline std::mt19937 benchNdGen{42};
inline std::uniform_real_distribution<double> benchNdDist(-1.0, 1.0);

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

void runBenchNdBuild(ankerl::nanobench::Bench &bench);
void runBenchNdEval(ankerl::nanobench::Bench &bench);
void runBenchNdMisc(ankerl::nanobench::Bench &bench);

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalSingle(const std::string &label, ankerl::nanobench::Bench &bench) {
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;

    Input x{};
    for (auto &xi : x) xi = benchNdDist(benchNdGen);

    const auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, unitCubeMin<Input>(), unitCubeMax<Input>());

    bench.run(label, [&] { ankerl::nanobench::doNotOptimizeAway(approx(x)); });
}

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS>
void benchEvalGenericPoint(const std::string &label, ankerl::nanobench::Bench &bench) {
    using CanonicalInput = std::array<double, IN_DIM>;
    using Input = FixedVec<double, IN_DIM>;
    using Output = FixedVec<double, OUT_DIM>;

    Input x{};
    for (std::size_t i = 0; i < IN_DIM; ++i) x[i] = benchNdDist(benchNdGen);

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
    for (auto &xi : x) xi = benchNdDist(benchNdGen);

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
        for (auto &xi : pt) xi = benchNdDist(benchNdGen);
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
        for (std::size_t i = 0; i < IN_DIM; ++i) pt[i] = benchNdDist(benchNdGen);
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

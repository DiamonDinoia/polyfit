#include <cmath>
#include <nanobench.h>
#include <random>
#include <vector>

#include "polyfit/polyfit.hpp"

// Helper to build a FuncEvalMany with N identical sin evaluators
template<std::size_t N, std::size_t... Is> static auto makeGroupImpl(std::index_sequence<Is...>) {
    auto makeOne = [] {
        return poly_eval::make_func_eval([](double x) { return std::sin(x); }, 16, -1.0, 1.0);
    };
    return poly_eval::make_func_eval_many((static_cast<void>(Is), makeOne())...);
}

template<std::size_t N> static auto makeGroup() { return makeGroupImpl<N>(std::make_index_sequence<N>{}); }

// Benchmark the “many” (tuple-packed) version
template<std::size_t N> static void benchGroup(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto group = makeGroup<N>();
    bench.run(std::to_string(N) + " funcs", [&] {
        for (double x : pts) {
            auto tup = group(x);
            ankerl::nanobench::doNotOptimizeAway(tup);
        }
    });
}

// Helper to build an array of N independent sin evaluators
template<std::size_t N, std::size_t... Is> static auto makeFuncsImpl(std::index_sequence<Is...>) {
    auto makeOne = [] {
        return poly_eval::make_func_eval([](double x) { return std::sin(x); }, 16, -1.0, 1.0);
    };
    return std::array{(static_cast<void>(Is), makeOne())...};
}

template<std::size_t N> static auto makeFuncs() { return makeFuncsImpl<N>(std::make_index_sequence<N>{}); }

// Benchmark the “non-many” version: loop over each func and eval x
template<std::size_t N> static void benchNonMany(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto funcs = makeFuncs<N>();
    bench.run(std::to_string(N) + " funcs (non-many)", [&] {
        for (double x : pts) {
            for (auto &f : funcs) {
                auto y = f(x);
                ankerl::nanobench::doNotOptimizeAway(y);
            }
        }
    });
}

// Compile-time dispatcher up to MaxNCoeffs, invoking both benchmarks
template<std::size_t MaxNCoeffs>
static bool dispatchBySize(std::size_t nCoeffs, const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    if constexpr (MaxNCoeffs == 0) {
        return false;
    } else {
        if (nCoeffs == MaxNCoeffs) {
            // tuple-packed “many” version
            benchGroup<MaxNCoeffs>(pts, bench);
            benchNonMany<MaxNCoeffs>(pts, bench);
            return true;
        }
        return dispatchBySize<MaxNCoeffs - 1>(nCoeffs, pts, bench);
    }
}

int main(int /*argc*/, char ** /*argv*/) {
    // 1. generate the inputs once
    constexpr std::size_t num_points = 1024;
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::vector<double> pts(num_points);
    for (auto &p : pts) p = dist(rng);

    // 2. configure the benchmark object once
    ankerl::nanobench::Bench bench;
    bench.title("poly_eval grouped-function throughput")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(10))
        .batch(num_points);

    // 3. run benchmarks for N = 1 .. 16
    for (std::size_t nCoeffs = 1; nCoeffs <= 16; ++nCoeffs) {
        dispatchBySize<16>(nCoeffs, pts, bench); // always succeeds for 1–16
    }

    return 0;
}

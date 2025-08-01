#include <cmath>
#include <nanobench.h>
#include <random>
#include <vector>

#include "polyfit/fast_eval.hpp"

// Helper to build a FuncEvalMany with N identical sin evaluators
template <std::size_t N, std::size_t... Is> static auto make_group_impl(std::index_sequence<Is...>) {
    auto make_one = [] { return poly_eval::make_func_eval([](double x) { return std::sin(x); }, 16, -1.0, 1.0); };
    return poly_eval::make_func_eval_many((static_cast<void>(Is), make_one())...);
}

template <std::size_t N> static auto make_group() { return make_group_impl<N>(std::make_index_sequence<N>{}); }

// Benchmark the “many” (tuple-packed) version
template <std::size_t N> static void bench_group(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto group = make_group<N>();
    bench.run(std::to_string(N) + " funcs", [&] {
        for (double x : pts) {
            auto tup = group(x);
            ankerl::nanobench::doNotOptimizeAway(tup);
        }
    });
}

// Helper to build an array of N independent sin evaluators
template <std::size_t N, std::size_t... Is> static auto make_funcs_impl(std::index_sequence<Is...>) {
    auto make_one = [] { return poly_eval::make_func_eval([](double x) { return std::sin(x); }, 16, -1.0, 1.0); };
    // build std::array of N copies
    return std::array{(static_cast<void>(Is), make_one())...};
}

template <std::size_t N> static auto make_funcs() { return make_funcs_impl<N>(std::make_index_sequence<N>{}); }

// Benchmark the “non-many” version: loop over each func and eval x
template <std::size_t N> static void bench_non_many(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto funcs = make_funcs<N>();
    bench.run(std::to_string(N) + " funcs (non-many)", [&] {
        for (double x : pts) {
            for (auto &f : funcs) {
                auto y = f(x);
                ankerl::nanobench::doNotOptimizeAway(y);
            }
        }
    });
}

// Compile-time dispatcher up to MaxN, invoking both benchmarks
template <std::size_t MaxN>
static bool dispatch(std::size_t n, const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    if constexpr (MaxN == 0) {
        return false;
    } else {
        if (n == MaxN) {
            // tuple-packed “many” version
            bench_group<MaxN>(pts, bench);
            // simple array-iteration version
            bench_non_many<MaxN>(pts, bench);
            return true;
        }
        return dispatch<MaxN - 1>(n, pts, bench);
    }
}

int main(int argc, char **argv) {
    // 1. generate the inputs once
    constexpr std::size_t num_points = 1024;
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::vector<double> pts(num_points);
    for (auto &p : pts)
        p = dist(rng);

    // 2. configure the benchmark object once
    ankerl::nanobench::Bench bench;
    bench.title("poly_eval grouped-function throughput")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(10))
        .batch(num_points);

    // 3. run benchmarks for N = 1 .. 16
    for (std::size_t n = 1; n <= 16; ++n) {
        dispatch<16>(n, pts, bench); // always succeeds for 1–16
    }

    return 0;
}

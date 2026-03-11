#include <cmath>
#include <nanobench.h>
#include <random>
#include <vector>

#include "polyfit/polyfit.hpp"

template<std::size_t N, class Builder, std::size_t... Is>
static auto makeArrayImpl(Builder build, std::index_sequence<Is...>) {
    return std::array{(static_cast<void>(Is), build())...};
}

template<std::size_t N> static auto makeEvalArray() {
    auto makeOne = [] {
        return poly_eval::fit([](double x) { return std::sin(x); }, 16, -1.0, 1.0);
    };
    return makeArrayImpl<N>(makeOne, std::make_index_sequence<N>{});
}

template<std::size_t N> static auto makeGroup() {
    return std::apply([](const auto &...evals) { return poly_eval::pack(evals...); }, makeEvalArray<N>());
}

template<std::size_t N> static void benchGroup(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto group = makeGroup<N>();
    bench.run(std::to_string(N) + " funcs", [&] {
        for (double x : pts) {
            auto tup = group(x);
            ankerl::nanobench::doNotOptimizeAway(tup);
        }
    });
}

template<std::size_t N> static void benchNonMany(const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    auto funcs = makeEvalArray<N>();
    bench.run(std::to_string(N) + " funcs (non-many)", [&] {
        for (double x : pts) {
            for (auto &f : funcs) {
                auto y = f(x);
                ankerl::nanobench::doNotOptimizeAway(y);
            }
        }
    });
}

template<std::size_t MaxNCoeffs>
static bool dispatchBySize(std::size_t nCoeffs, const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    if constexpr (MaxNCoeffs == 0) {
        return false;
    } else {
        if (nCoeffs == MaxNCoeffs) {
            benchGroup<MaxNCoeffs>(pts, bench);
            benchNonMany<MaxNCoeffs>(pts, bench);
            return true;
        }
        return dispatchBySize<MaxNCoeffs - 1>(nCoeffs, pts, bench);
    }
}

int main(int /*argc*/, char ** /*argv*/) {
    constexpr std::size_t numPoints = 1024;
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::vector<double> pts(numPoints);
    for (auto &p : pts) p = dist(rng);

    ankerl::nanobench::Bench bench;
    bench.title("poly_eval grouped-function throughput")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(10))
        .batch(numPoints);

    for (std::size_t nCoeffs = 1; nCoeffs <= 16; ++nCoeffs) {
        dispatchBySize<16>(nCoeffs, pts, bench);
    }

    return 0;
}

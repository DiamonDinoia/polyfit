#include <array>
#include <cmath>
#include <complex>
#include <nanobench.h>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "polyfit/polyfit.hpp"

static double bench_func(double x) { return std::sin(x); }
static std::complex<double> bench_complex_func(double x) { return {std::sin(x), std::cos(x)}; }

template<std::size_t I> static auto makeMixedEval() {
    if constexpr ((I % 2) == 0) {
        return poly_eval::fit(bench_func, 16, -1.0, 1.0);
    } else {
        return poly_eval::fit(bench_func, 16, 0.0, 1.0);
    }
}

template<std::size_t N, class Builder, std::size_t... Is>
static auto makeArrayImpl(Builder build, std::index_sequence<Is...>) {
    return std::array{(static_cast<void>(Is), build())...};
}

template<std::size_t N> static auto makeEvalArray() {
    auto makeOne = [] {
        return poly_eval::fit(bench_func, 16, -1.0, 1.0);
    };
    return makeArrayImpl<N>(makeOne, std::make_index_sequence<N>{});
}

template<std::size_t N, std::size_t... Is> static auto makeMixedEvalArrayImpl(std::index_sequence<Is...>) {
    return std::array{makeMixedEval<Is>()...};
}

template<std::size_t N> static auto makeMixedEvalArray() {
    return makeMixedEvalArrayImpl<N>(std::make_index_sequence<N>{});
}

template<std::size_t N> static auto makeGroup() {
    return std::apply([](const auto &...evals) { return poly_eval::pack(evals...); }, makeEvalArray<N>());
}

template<std::size_t N> static auto makeMixedGroup() {
    return std::apply([](const auto &...evals) { return poly_eval::pack(evals...); }, makeMixedEvalArray<N>());
}

template<std::size_t N> static auto makeComplexEvalArray() {
    auto makeOne = [] {
        return poly_eval::fit(bench_complex_func, 16, -1.0, 1.0);
    };
    return makeArrayImpl<N>(makeOne, std::make_index_sequence<N>{});
}

template<std::size_t N> static auto makeComplexGroup() {
    return std::apply([](const auto &...evals) { return poly_eval::pack(evals...); }, makeComplexEvalArray<N>());
}

template<std::size_t N, class Group>
static void benchGroup(const std::vector<double> &pts, ankerl::nanobench::Bench &bench, Group &&group,
                       const std::string &name) {
    bench.run(name, [&] {
        for (double x : pts) {
            auto tup = group(x);
            ankerl::nanobench::doNotOptimizeAway(tup);
        }
    });
}

template<std::size_t N, class Evals>
static void benchNonMany(const std::vector<double> &pts, ankerl::nanobench::Bench &bench, Evals &&funcs,
                         const std::string &name) {
    bench.run(name, [&] {
        for (double x : pts) {
            for (auto &f : funcs) {
                auto y = f(x);
                ankerl::nanobench::doNotOptimizeAway(y);
            }
        }
    });
}

template<std::size_t N, class Group>
static void benchGroupBulk(const std::vector<double> &pts, ankerl::nanobench::Bench &bench, Group &&group,
                           const std::string &name) {
    std::vector<typename std::remove_reference_t<Group>::OutputType> out(pts.size() * N);
    bench.run(name, [&] {
        group(pts.data(), out.data(), pts.size());
        ankerl::nanobench::doNotOptimizeAway(out.data());
    });
}

template<std::size_t MaxNCoeffs>
static bool dispatchBySize(std::size_t nCoeffs, const std::vector<double> &pts, ankerl::nanobench::Bench &bench) {
    if constexpr (MaxNCoeffs == 0) {
        return false;
    } else {
        if (nCoeffs == MaxNCoeffs) {
            benchGroup<MaxNCoeffs>(pts, bench, makeGroup<MaxNCoeffs>(), std::to_string(MaxNCoeffs) + " funcs");
            benchNonMany<MaxNCoeffs>(pts, bench, makeEvalArray<MaxNCoeffs>(),
                                     std::to_string(MaxNCoeffs) + " funcs (non-many)");
            return true;
        }
        return dispatchBySize<MaxNCoeffs - 1>(nCoeffs, pts, bench);
    }
}

int main(int /*argc*/, char ** /*argv*/) {
    constexpr std::size_t numPoints = 1024;
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::uniform_real_distribution<double> mixed_dist(0.0, 1.0);

    std::vector<double> pts(numPoints);
    for (auto &p : pts) p = dist(rng);
    std::vector<double> mixed_pts(numPoints);
    for (auto &p : mixed_pts) p = mixed_dist(rng);

    ankerl::nanobench::Bench bench;
    bench.title("poly_eval grouped-function throughput")
        .unit("eval")
        .minEpochTime(std::chrono::milliseconds(20))
        .minEpochIterations(100)
        .batch(numPoints);

    for (std::size_t nCoeffs = 1; nCoeffs <= 16; ++nCoeffs) {
        dispatchBySize<16>(nCoeffs, pts, bench);
    }

    benchGroup<8>(mixed_pts, bench, makeMixedGroup<8>(), "8 funcs (mixed domains)");
    benchNonMany<8>(mixed_pts, bench, makeMixedEvalArray<8>(), "8 funcs (mixed domains, non-many)");
    benchGroup<16>(mixed_pts, bench, makeMixedGroup<16>(), "16 funcs (mixed domains)");
    benchNonMany<16>(mixed_pts, bench, makeMixedEvalArray<16>(), "16 funcs (mixed domains, non-many)");
    benchGroupBulk<8>(pts, bench, makeComplexGroup<8>(), "8 funcs (complex outputs, bulk)");

    return 0;
}

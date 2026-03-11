#include <chrono>
#include <complex>
#include <nanobench.h>
#include <random>
#include <vector>

#include "polyfit/polyfit.hpp"

using namespace ankerl::nanobench;

std::mt19937_64 rng(42);

template<typename F>
void runBench(const std::string &label, Bench &bench, F func, typename poly_eval::FunctionTraits<F>::arg0_type a,
              typename poly_eval::FunctionTraits<F>::arg0_type b, size_t numPoints = 1024) {
    using T = decltype(a);
    using V = typename poly_eval::FunctionTraits<F>::result_type;

    std::uniform_real_distribution<T> dist(a, b);
    std::vector<T> pts;
    pts.reserve(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        if constexpr (std::is_floating_point_v<T>) {
            pts.push_back(dist(rng));
        } else {
            pts.emplace_back(dist(rng), dist(rng));
        }
    }
    std::vector<V> out(numPoints);

    bench.run(label + " constexpr fit",
              [&] { doNotOptimizeAway(poly_eval::fit<8>(func, a, b, poly_eval::Iters<0>{})); });

    bench.run(label + " coefficient-count fit", [&] { doNotOptimizeAway(poly_eval::fit(func, 8, a, b)); });

    bench.run(label + " eps fit", [&] { doNotOptimizeAway(poly_eval::fit(func, 1e-6, a, b)); });

    const auto fe = poly_eval::fit<8>(func, a, b, poly_eval::Iters<0>{});
    bench.run(label + " eval", [&] {
        auto r = fe(pts[0]);
        doNotOptimizeAway(r);
    });

    bench.batch(pts.size());
    bench.run(label + " eval_many", [&] {
        fe(pts.data(), out.data(), pts.size());
        doNotOptimizeAway(out);
    });
    bench.batch(1);
}

int main() {
    Bench bench;
    bench.minEpochTime(std::chrono::milliseconds(10));
    runBench("float", bench, [](float x) { return x * x + 1.0f; }, -1.5f, 1.5f);
    runBench("double", bench, [](double x) { return x * x + 1.0; }, -2.0, 2.0);
    runBench("complex<double>", bench, [](double x) { return x * x + std::complex<double>(1.0, 1.0); }, -1.5, 1.5);
    runBench("complex<float>", bench, [](float x) { return x * x + std::complex<float>(1.0, 1.0); }, -1.2f, 1.3f);
    return 0;
}

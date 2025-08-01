#include <chrono>
#include <complex>
#include <nanobench.h>
#include <random>
#include <vector>

#include "polyfit/fast_eval.hpp"

using namespace ankerl::nanobench;

std::mt19937_64 rng(42);
// A single templated bench function: T is the input/output type, F is the functor.
template <typename F>
void run_bench(const std::string &label, Bench &bench, F func, typename poly_eval::function_traits<F>::arg0_type a,
               typename poly_eval::function_traits<F>::arg0_type b, size_t num_points = 1024) {
    using T = decltype(a);
    using V = typename poly_eval::function_traits<F>::result_type;

    // Prepare random points and output buffer
    std::uniform_real_distribution<T> dist(a, b);
    std::vector<T> pts;
    pts.reserve(num_points);
    for (size_t i = 0; i < num_points; ++i) {
        if constexpr (std::is_floating_point_v<T>) {
            pts.push_back(dist(rng));
        } else {
            pts.emplace_back(dist(rng), dist(rng));
        }
    }
    std::vector<V> out(num_points);

    bench.run(label + " constexpr fit", [&] { doNotOptimizeAway(poly_eval::make_func_eval<8, 0>(func, a, b)); });

    bench.run(label + " degree fit", [&] { doNotOptimizeAway(poly_eval::make_func_eval(func, 8, a, b)); });

    bench.run(label + " eps fit", [&] { doNotOptimizeAway(poly_eval::make_func_eval(func, 1e-6, a, b)); });

    const auto fe = poly_eval::make_func_eval<8, 0>(func, a, b);
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
    run_bench("float", bench, [](float x) { return x * x + 1.0f; }, -1.5f, 1.5f);
    run_bench("double", bench, [](double x) { return x * x + 1.0; }, -2.0, 2.0);
    run_bench("complex<double>", bench, [](double x) { return x * x + std::complex<double>(1.0, 1.0); }, -1.5, 1.5);
    run_bench("complex<float>", bench, [](float x) { return x * x + std::complex<float>(1.0, 1.0); }, -1.2f, 1.3f);
    return 0;
}

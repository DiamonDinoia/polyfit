#include "polyfit/fast_eval.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>

int main() {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    const std::vector<int> unrolls = {1, 2, 4, 8, 16, 32};
    const std::vector<int> degs = {8, 16, 24, 32};

    const std::size_t P = 1024;            // points per eval
    const int warmup_iters = 10;
    const int bench_iters = 200;           // iterations to average over

    // storage for results: rows = unrolls, cols = degs
    std::vector<std::vector<double>> results(unrolls.size(), std::vector<double>(degs.size(), 0.0));

    std::vector<double> pts(P), out(P);

    for (auto &p : pts) p = dist(rng);

    std::cout << "| Unroll \\ Deg |";
    for (auto d : degs) std::cout << " " << d << " |";
    std::cout << "\n|---|";
    for (size_t i = 0; i < degs.size(); ++i) std::cout << "---|";
    std::cout << "\n";

    for (size_t ui = 0; ui < unrolls.size(); ++ui) {
        int U = unrolls[ui];
        std::cout << "| " << U << " |";
        for (size_t di = 0; di < degs.size(); ++di) {
            int deg = degs[di];

            // prepare coefficients
            std::vector<double> coeffs(deg);
            for (auto &c : coeffs) c = dist(rng);

            // warmup
            for (int w = 0; w < warmup_iters; ++w) {
                poly_eval::horner<0, false, false, 0>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                     [](auto v) { return v; }); // baseline call to warm caches
            }
            // timed runs: use the templated UNROLL parameter by switching on U
            auto run_templated = [&](int unroll) {
                using namespace std::chrono;
                auto t0 = steady_clock::now();
                for (int it = 0; it < bench_iters; ++it) {
                    switch (unroll) {
                    case 1:
                        poly_eval::horner<0, false, false, 1>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                             [](auto v) { return v; });
                        break;
                    case 2:
                        poly_eval::horner<0, false, false, 2>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                             [](auto v) { return v; });
                        break;
                    case 4:
                        poly_eval::horner<0, false, false, 4>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                             [](auto v) { return v; });
                        break;
                    case 8:
                        poly_eval::horner<0, false, false, 8>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                             [](auto v) { return v; });
                        break;
                    case 16:
                        poly_eval::horner<0, false, false, 16>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                              [](auto v) { return v; });
                        break;
                    case 32:
                        poly_eval::horner<0, false, false, 32>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                              [](auto v) { return v; });
                        break;
                    default:
                        // Fallback: call with 0 (auto) if unexpected
                        poly_eval::horner<0, false, false, 0>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                             [](auto v) { return v; });
                    }
                }
                auto t1 = steady_clock::now();
                return duration_cast<nanoseconds>(t1 - t0).count();
            };

            // measure
            auto ns_total = run_templated(U);
            double ns_per_eval = double(ns_total) / double(bench_iters * P);
            results[ui][di] = ns_per_eval;

            std::cout << " " << std::fixed << std::setprecision(2) << ns_per_eval << " |";
            // flush to keep progress visible in CI
            std::cout << std::flush;
        }
        std::cout << "\n";
    }

    // Print a compact compatibility report below the table
    std::cout << "\n## Unroll sweep summary (ns/eval)\n\n";
    std::cout << "| Unroll |";
    for (auto d : degs) std::cout << " " << d << " |";
    std::cout << "\n|---|";
    for (size_t i = 0; i < degs.size(); ++i) std::cout << "---|";
    std::cout << "\n";
    for (size_t ui = 0; ui < unrolls.size(); ++ui) {
        std::cout << "| " << unrolls[ui] << " |";
        for (size_t di = 0; di < degs.size(); ++di) {
            std::cout << " " << std::fixed << std::setprecision(2) << results[ui][di] << " |";
        }
        std::cout << "\n";
    }

    return 0;
}

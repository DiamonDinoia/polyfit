// bench_horner.cpp
// -----------------------------------------------------------------------------
// NanoBench suite for poly_eval public APIs with per‑point timing (dynamic only):
// This replaces the Google Benchmark driver entirely—ensure CMakeLists points to
// this file
//   - Scalar Horner (runtime)
//   - SIMD Horner (runtime)
//   - horner_many (runtime)
//   - ND Horner front‑end (2D, 3D, 4D)
//
// Uses .batch(Points) and .context("FMAs", ...) to label FMA counts, then renders
// a Markdown table that automatically includes an `FMAs` column via
// ankerl::nanobench::templates::markdown().
// -----------------------------------------------------------------------------

 
#include "polyfit/fast_eval.hpp"
#include <nanobench.h>

#include <array>
#include <cmath>
#include <experimental/mdspan>
#include <random>
#include <string>
#include <vector>

int main() {
    // RNG setup
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Benchmark configuration
    ankerl::nanobench::Bench bench;
    bench.title("poly_eval Horner Suite");
    bench.unit("eval");
    bench.epochs(50);

    // 1D: scalar and SIMD runtime for degrees 8,16,24,32
    std::vector<int> degs1d = {8, 16, 24, 32};
    for (auto deg : degs1d) {
        // Scalar Horner runtime
        bench.context("FMAs", std::to_string(deg))
            .batch(1)
            .run("Dim=1, Deg=" + std::to_string(deg) + ", SIMD=No, scalar-runtime",
                 [&] {
                     std::vector<double> coeffs(deg);
                     for (auto &c : coeffs)
                         c = dist(rng);
                     double x = dist(rng);
                     ankerl::nanobench::doNotOptimizeAway(poly_eval::horner(x, coeffs.data(), deg));
                 })
            .clearContext();

        // SIMD Horner runtime
        bench.context("FMAs", std::to_string(deg))
            .batch(1024)
            .run("Dim=1, Deg=" + std::to_string(deg) + ", SIMD=Yes, simd-runtime",
                 [&] {
                     const size_t P = 1024;
                     std::vector<double> pts(P), out(P), coeffs(deg);
                     for (auto &p : pts)
                         p = dist(rng);
                     for (auto &c : coeffs)
                         c = dist(rng);
                     poly_eval::horner<0, false, false, 0>(pts.data(), out.data(), P, coeffs.data(), deg,
                                                           [](auto v) { return v; });
                     ankerl::nanobench::doNotOptimizeAway(out.data());
                 })
            .clearContext();
    }

    // horner_many runtime for various (M, Deg)
    std::vector<int> Ms = {4, 8, 12, 16};
    for (auto deg : degs1d) {
        for (auto M : Ms) {
            bench.context("FMAs", std::to_string(M * deg))
                .batch(M)
                .run("Dim=1, horner_many M=" + std::to_string(M) + ", Deg=" + std::to_string(deg) +
                         ", SIMD=No, many-runtime",
                     [&] {
                         std::vector<double> coeffs(M * deg), out(M);
                         for (auto &c : coeffs)
                             c = dist(rng);
                         double x = dist(rng);
                         poly_eval::horner_many<0, 0, false, double, double>(x, coeffs.data(), out.data(), M, deg);
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();
        }
    }

    // 2D ND Horner runtime for degrees 4,8,16,24
    std::vector<int> degs2d = {4, 8, 16, 24};
    for (auto deg : degs2d) {
        bench.context("FMAs", std::to_string(deg * deg * 2))
            .batch(1)
            .run("Dim=2, Deg=" + std::to_string(deg) + ", SIMD=No, ND2D",
                 [&] {
                     const size_t Dim = 2;
                     std::vector<double> coeffs(deg * deg * Dim);
                     for (auto &c : coeffs)
                         c = dist(rng);
                     std::experimental::mdspan<const double, std::experimental::dextents<size_t, 3>> C(coeffs.data(),
                                                                                                       deg, deg, Dim);
                     std::array<double, Dim> x{dist(rng), dist(rng)};
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(x, C, deg));
                 })
            .clearContext();
    }

    // 3D ND Horner runtime for degrees 4,8,16
    std::vector<int> degs3d = {4, 8, 16};
    for (auto deg : degs3d) {
        bench.context("FMAs", std::to_string(deg * deg * deg * 3))
            .batch(1)
            .run("Dim=3, Deg=" + std::to_string(deg) + ", SIMD=No, ND3D",
                 [&] {
                     const size_t Dim = 3;
                     std::vector<double> coeffs(deg * deg * deg * Dim);
                     for (auto &c : coeffs)
                         c = dist(rng);
                     std::experimental::mdspan<const double, std::experimental::dextents<size_t, 4>> C(
                         coeffs.data(), deg, deg, deg, Dim);
                     std::array<double, Dim> x{dist(rng), dist(rng), dist(rng)};
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(x, C, deg));
                 })
            .clearContext();
    }

    // 4D ND Horner runtime for degrees 4,8
    std::vector<int> degs4d = {4, 8};
    for (auto deg : degs4d) {
        bench.context("FMAs", std::to_string(static_cast<int>(std::pow(deg, 4) * 4)))
            .batch(1)
            .run("Dim=4, Deg=" + std::to_string(deg) + ", SIMD=No, ND4D",
                 [&] {
                     const size_t Dim = 4;
                     size_t total = static_cast<size_t>(std::pow(deg, 4) * Dim);
                     std::vector<double> coeffs(total);
                     for (auto &c : coeffs)
                         c = dist(rng);
                     std::experimental::mdspan<const double, std::experimental::dextents<size_t, 5>> C(
                         coeffs.data(), deg, deg, deg, deg, Dim);
                     std::array<double, Dim> x{dist(rng), dist(rng), dist(rng), dist(rng)};
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(x, C, deg));
                 })
            .clearContext();
    }

    return 0;
}

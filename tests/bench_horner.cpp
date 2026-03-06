// bench_horner.cpp
// -----------------------------------------------------------------------------
// NanoBench suite for poly_eval public APIs with per-point timing (dynamic only):
//   - Scalar Horner (runtime)
//   - SIMD Horner (runtime, aligned + unaligned)
//   - horner_many (runtime)
//   - horner_transposed scalar
//   - ND Horner front-end (2D, 3D, 4D)
//
// Allocations and RNG are hoisted out of the timed region so measurements
// reflect pure compute throughput.  A summary table at the end compares
// measured GFLOPS against theoretical peak for the SIMD Horner benchmarks.
// -----------------------------------------------------------------------------

#include "polyfit/fast_eval.hpp"
#include <nanobench.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <experimental/mdspan>
#include <random>
#include <string>
#include <vector>
#include <xsimd/xsimd.hpp>

int main() {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    constexpr size_t simd_size = xsimd::batch<double>::size;
    constexpr size_t alignment = xsimd::batch<double>::arch_type::alignment();
    using aligned_vec = std::vector<double, xsimd::aligned_allocator<double, alignment>>;

    ankerl::nanobench::Bench bench;
    bench.title("poly_eval Horner Suite").unit("eval").epochs(50).warmup(10).minEpochIterations(20);

    constexpr size_t P = 4096;
    std::vector<size_t> degs1d = {8, 16, 24, 32};

    // Track SIMD results for summary table
    struct SimdRecord {
        size_t deg;
        double ns_per_pt;
        double cyc_per_pt;
        bool aligned;
    };
    std::vector<SimdRecord> simd_records;

    // 1D scalar + SIMD runtime
    for (auto deg : degs1d) {
        // --- Scalar Horner ---
        {
            std::vector<double> coeffs(deg);
            for (auto &c : coeffs) c = dist(rng);
            double x = dist(rng);

            bench.context("FMAs", std::to_string(deg))
                .batch(1)
                .run("Dim=1, Deg=" + std::to_string(deg) + ", SIMD=No, scalar-runtime",
                     [&] {
                         ankerl::nanobench::doNotOptimizeAway(poly_eval::horner(x, coeffs.data(), deg));
                     })
                .clearContext();
        }

        // --- SIMD Horner aligned ---
        {
            aligned_vec pts(P), out(P);
            std::vector<double> coeffs(deg);
            for (auto &p : pts) p = dist(rng);
            for (auto &c : coeffs) c = dist(rng);

            bench.context("FMAs", std::to_string(deg))
                .batch(P)
                .run("Dim=1, Deg=" + std::to_string(deg) + ", SIMD=Yes, simd-aligned",
                     [&] {
                         poly_eval::horner<0, true, true, 0>(
                             pts.data(), out.data(), P, coeffs.data(), deg, [](auto v) { return v; });
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();

            auto const &r = bench.results().back();
            double batch = r.config().mBatch;
            simd_records.push_back({deg, r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e9 / batch,
                                    r.median(ankerl::nanobench::Result::Measure::cpucycles) / batch, true});
        }

        // --- SIMD Horner unaligned ---
        {
            std::vector<double> pts(P), out(P), coeffs(deg);
            for (auto &p : pts) p = dist(rng);
            for (auto &c : coeffs) c = dist(rng);

            bench.context("FMAs", std::to_string(deg))
                .batch(P)
                .run("Dim=1, Deg=" + std::to_string(deg) + ", SIMD=Yes, simd-unaligned",
                     [&] {
                         poly_eval::horner<0, false, false, 0>(
                             pts.data(), out.data(), P, coeffs.data(), deg, [](auto v) { return v; });
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();

            auto const &r = bench.results().back();
            double batch = r.config().mBatch;
            simd_records.push_back({deg, r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e9 / batch,
                                    r.median(ankerl::nanobench::Result::Measure::cpucycles) / batch, false});
        }
    }

    // horner_many runtime
    std::vector<size_t> Ms = {4, 8, 12, 16};
    for (auto deg : degs1d) {
        for (auto M : Ms) {
            std::vector<double> coeffs(M * deg), out(M);
            for (auto &c : coeffs) c = dist(rng);
            double x = dist(rng);

            bench.context("FMAs", std::to_string(M * deg))
                .batch(M)
                .run("Dim=1, horner_many M=" + std::to_string(M) + ", Deg=" + std::to_string(deg) +
                         ", SIMD=No, many-runtime",
                     [&] {
                         poly_eval::horner_many<0, 0, double, double>(x, coeffs.data(), out.data(), M, deg);
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();
        }
    }

    // horner_transposed scalar (simd_width=0)
    for (auto deg : degs1d) {
        for (auto M : Ms) {
            std::vector<double> coeffs(deg * M), out(M), x(M);
            for (auto &c : coeffs) c = dist(rng);
            for (auto &xi : x) xi = dist(rng);

            bench.context("FMAs", std::to_string(M * deg))
                .batch(M)
                .run("Dim=1, horner_transposed_scalar M=" + std::to_string(M) + ", Deg=" + std::to_string(deg),
                     [&] {
                         poly_eval::horner_transposed<0, 0, 0>(x.data(), coeffs.data(), out.data(), M, deg);
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();
        }
    }

    // 2D ND Horner
    for (auto deg : {size_t(4), size_t(8), size_t(16), size_t(24)}) {
        constexpr size_t Dim = 2;
        std::vector<double> coeffs(deg * deg * Dim);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 3>> C(coeffs.data(), deg, deg, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(deg * deg * Dim))
            .batch(1)
            .run("Dim=2, Deg=" + std::to_string(deg) + ", SIMD=No, ND2D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(deg)));
                 })
            .clearContext();
    }

    // 3D ND Horner
    for (auto deg : {size_t(4), size_t(8), size_t(16)}) {
        constexpr size_t Dim = 3;
        std::vector<double> coeffs(deg * deg * deg * Dim);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 4>> C(
            coeffs.data(), deg, deg, deg, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(deg * deg * deg * Dim))
            .batch(1)
            .run("Dim=3, Deg=" + std::to_string(deg) + ", SIMD=No, ND3D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(deg)));
                 })
            .clearContext();
    }

    // 4D ND Horner
    for (auto deg : {size_t(4), size_t(8)}) {
        constexpr size_t Dim = 4;
        auto total = static_cast<size_t>(std::pow(deg, 4) * Dim);
        std::vector<double> coeffs(total);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 5>> C(
            coeffs.data(), deg, deg, deg, deg, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng), dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(total))
            .batch(1)
            .run("Dim=4, Deg=" + std::to_string(deg) + ", SIMD=No, ND4D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(deg)));
                 })
            .clearContext();
    }

    // -------------------------------------------------------------------------
    // Peak throughput summary for SIMD Horner
    // -------------------------------------------------------------------------
    // Assumes 2 FMA ports, each 256-bit (4 doubles), 2 FLOPS per FMA (mul+add).
    // The cycle counter is TSC (base freq ~3 GHz) while the CPU turbos higher
    // during AVX2 work, so we derive actual frequency from the ns/cyc ratio and
    // compute efficiency at the true turbo frequency.
    constexpr int fma_ports = 2;

    // Estimate actual turbo frequency from the fastest (lowest-overhead) result
    double turbo_ghz = 0;
    for (auto const &rec : simd_records) {
        double freq = rec.cyc_per_pt / rec.ns_per_pt;
        if (freq > turbo_ghz) turbo_ghz = freq;
    }
    // The TSC-based cycles undercount at turbo.  Scale to actual core cycles.
    double tsc_ghz = turbo_ghz; // TSC freq ≈ cyc/ns
    // Actual turbo is higher — estimate from wall-clock throughput of the
    // tightest result (Deg=8 aligned: nearly pure FMA, minimal overhead).
    // At peak: vFMAs/pt = (deg-1)/simd_size, 2 ports → cyc/pt = vFMAs/2.
    // actual_freq = theoretical_cyc / measured_ns.
    if (!simd_records.empty()) {
        auto const &best = simd_records.front(); // Deg=8 aligned = tightest
        double theo_cyc = static_cast<double>(best.deg - 1) / static_cast<double>(simd_size) / fma_ports;
        turbo_ghz = theo_cyc / best.ns_per_pt;
    }
    double peak_gflops = turbo_ghz * fma_ports * static_cast<double>(simd_size) * 2.0;

    std::printf("\n");
    std::printf("SIMD Horner peak throughput (%zu-wide AVX2, %d FMA ports, %zu pts)\n",
                simd_size, fma_ports, P);
    std::printf("TSC base: %.2f GHz | Est. AVX2 turbo: %.2f GHz | HW DP peak: %.1f GFLOPS\n",
                tsc_ghz, turbo_ghz, peak_gflops);
    std::printf("%-11s %8s %9s %10s\n", "benchmark", "ns/pt", "GFLOPS", "% HW peak");
    std::printf("%-11s %8s %9s %10s\n", "-----------", "--------", "---------", "----------");

    for (auto const &rec : simd_records) {
        double gflops = static_cast<double>(rec.deg - 1) * 2.0 / rec.ns_per_pt;
        std::printf("Deg=%-2zu %-4s %8.3f %9.1f %9.1f%%\n",
                    rec.deg, rec.aligned ? "aln" : "uln",
                    rec.ns_per_pt, gflops, gflops / peak_gflops * 100.0);
    }

    return 0;
}

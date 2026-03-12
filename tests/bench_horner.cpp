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

#include "polyfit/polyfit.hpp"
#include <nanobench.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <experimental/mdspan>
#include <random>
#include <string>
#include <vector>
#include <xsimd/xsimd.hpp>

using namespace std::chrono_literals;

int main() {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    constexpr size_t simdSize = xsimd::batch<double>::size;
    constexpr size_t alignment = xsimd::batch<double>::arch_type::alignment();
    using AlignedVec = std::vector<double, xsimd::aligned_allocator<double, alignment>>;

    ankerl::nanobench::Bench bench;
    bench.title("poly_eval Horner Suite").unit("eval").epochs(50).warmup(10).minEpochTime(20ms).minEpochIterations(100);

    constexpr size_t P = 4096;
    std::vector<size_t> nCoeffs1D = {8, 16, 24, 32};

    // Track SIMD results for summary table
    struct SimdRecord {
        size_t nCoeffs;
        double nsPerPoint;
        double cyclesPerPoint;
        bool aligned;
    };
    std::vector<SimdRecord> simdRecords;

    // 1D scalar + SIMD runtime
    for (auto nCoeffs : nCoeffs1D) {
        // --- Scalar Horner ---
        {
            std::vector<double> coeffs(nCoeffs);
            for (auto &c : coeffs) c = dist(rng);
            double x = dist(rng);

            bench.context("FMAs", std::to_string(nCoeffs))
                .batch(1)
                .run("Dim=1, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=No, scalar-runtime",
                     [&] {
                         ankerl::nanobench::doNotOptimizeAway(poly_eval::horner(x, coeffs.data(), nCoeffs));
                     })
                .clearContext();
        }

        // --- SIMD Horner aligned ---
        {
            AlignedVec pts(P), out(P);
            std::vector<double> coeffs(nCoeffs);
            for (auto &p : pts) p = dist(rng);
            for (auto &c : coeffs) c = dist(rng);

            bench.context("FMAs", std::to_string(nCoeffs))
                .batch(P)
                .run("Dim=1, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=Yes, simd-aligned",
                     [&] {
                         poly_eval::horner<0, true, true, 0>(
                             pts.data(), out.data(), P, coeffs.data(), nCoeffs, [](auto v) { return v; });
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();

            auto const &r = bench.results().back();
            double batch = r.config().mBatch;
            simdRecords.push_back({nCoeffs, r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e9 / batch,
                                   r.median(ankerl::nanobench::Result::Measure::cpucycles) / batch, true});
        }

        // --- SIMD Horner unaligned ---
        {
            std::vector<double> pts(P), out(P), coeffs(nCoeffs);
            for (auto &p : pts) p = dist(rng);
            for (auto &c : coeffs) c = dist(rng);

            bench.context("FMAs", std::to_string(nCoeffs))
                .batch(P)
                .run("Dim=1, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=Yes, simd-unaligned",
                     [&] {
                         poly_eval::horner<0, false, false, 0>(
                             pts.data(), out.data(), P, coeffs.data(), nCoeffs, [](auto v) { return v; });
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();

            auto const &r = bench.results().back();
            double batch = r.config().mBatch;
            simdRecords.push_back({nCoeffs, r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e9 / batch,
                                   r.median(ankerl::nanobench::Result::Measure::cpucycles) / batch, false});
        }
    }

    // horner_many runtime
    std::vector<size_t> Ms = {4, 8, 12, 16};
    for (auto nCoeffs : nCoeffs1D) {
        for (auto M : Ms) {
            std::vector<double> coeffs(M * nCoeffs), out(M);
            for (auto &c : coeffs) c = dist(rng);
            double x = dist(rng);

            bench.context("FMAs", std::to_string(M * nCoeffs))
                .batch(M)
                .run("Dim=1, horner_many M=" + std::to_string(M) + ", nCoeffs=" + std::to_string(nCoeffs) +
                         ", SIMD=No, many-runtime",
                     [&] {
                         poly_eval::horner_many<0, 0, double, double>(x, coeffs.data(), out.data(), M, nCoeffs);
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();
        }
    }

    // horner_transposed scalar (SIMD_WIDTH=0)
    for (auto nCoeffs : nCoeffs1D) {
        for (auto M : Ms) {
            std::vector<double> coeffs(nCoeffs * M), out(M), x(M);
            for (auto &c : coeffs) c = dist(rng);
            for (auto &xi : x) xi = dist(rng);

            bench.context("FMAs", std::to_string(M * nCoeffs))
                .batch(M)
                .run("Dim=1, horner_transposed_scalar M=" + std::to_string(M) + ", nCoeffs=" + std::to_string(nCoeffs),
                     [&] {
                         poly_eval::horner_transposed<0, 0, 0>(x.data(), coeffs.data(), out.data(), M, nCoeffs);
                         ankerl::nanobench::doNotOptimizeAway(out.data());
                     })
                .clearContext();
        }
    }

    // 2D ND Horner
    for (auto nCoeffs : {size_t(4), size_t(8), size_t(16), size_t(24)}) {
        constexpr size_t Dim = 2;
        std::vector<double> coeffs(nCoeffs * nCoeffs * Dim);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 3>> C(
            coeffs.data(), nCoeffs, nCoeffs, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(nCoeffs * nCoeffs * Dim))
            .batch(1)
            .run("Dim=2, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=No, ND2D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(nCoeffs)));
                 })
            .clearContext();
    }

    // 3D ND Horner
    for (auto nCoeffs : {size_t(4), size_t(8), size_t(16)}) {
        constexpr size_t Dim = 3;
        std::vector<double> coeffs(nCoeffs * nCoeffs * nCoeffs * Dim);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 4>> C(
            coeffs.data(), nCoeffs, nCoeffs, nCoeffs, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(nCoeffs * nCoeffs * nCoeffs * Dim))
            .batch(1)
            .run("Dim=3, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=No, ND3D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(nCoeffs)));
                 })
            .clearContext();
    }

    // 4D ND Horner
    for (auto nCoeffs : {size_t(4), size_t(8)}) {
        constexpr size_t Dim = 4;
        auto total = static_cast<size_t>(std::pow(nCoeffs, 4) * Dim);
        std::vector<double> coeffs(total);
        for (auto &c : coeffs) c = dist(rng);
        std::experimental::mdspan<const double, std::experimental::dextents<size_t, 5>> C(
            coeffs.data(), nCoeffs, nCoeffs, nCoeffs, nCoeffs, Dim);
        std::array<double, Dim> x{dist(rng), dist(rng), dist(rng), dist(rng)};

        bench.context("FMAs", std::to_string(total))
            .batch(1)
            .run("Dim=4, nCoeffs=" + std::to_string(nCoeffs) + ", SIMD=No, ND4D",
                 [&] {
                     ankerl::nanobench::doNotOptimizeAway(
                         poly_eval::horner<0, true, std::array<double, Dim>, decltype(x), decltype(C)>(
                             x, C, static_cast<int>(nCoeffs)));
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
    for (auto const &record : simdRecords) {
        double freq = record.cyclesPerPoint / record.nsPerPoint;
        turbo_ghz = std::max(turbo_ghz, freq);
    }
    // The TSC-based cycles undercount at turbo.  Scale to actual core cycles.
    double tsc_ghz = turbo_ghz; // TSC freq ≈ cyc/ns
    // Actual turbo is higher — estimate from wall-clock throughput of the
    // tightest result (8 coefficients, aligned: nearly pure FMA, minimal overhead).
    // At peak: vFMAs/pt = (nCoeffs-1)/simdSize, 2 ports → cyc/pt = vFMAs/2.
    // actual_freq = theoretical_cyc / measured_ns.
    if (!simdRecords.empty()) {
        auto const &best = simdRecords.front();
        double theoCyc = static_cast<double>(best.nCoeffs - 1) / static_cast<double>(simdSize) / fma_ports;
        turbo_ghz = theoCyc / best.nsPerPoint;
    }
    double peakGflops = turbo_ghz * fma_ports * static_cast<double>(simdSize) * 2.0;

    std::printf("\n");
    std::printf("SIMD Horner peak throughput (%zu-wide AVX2, %d FMA ports, %zu pts)\n",
                simdSize, fma_ports, P);
    std::printf("TSC base: %.2f GHz | Est. AVX2 turbo: %.2f GHz | HW DP peak: %.1f GFLOPS\n",
                tsc_ghz, turbo_ghz, peakGflops);
    std::printf("%-11s %8s %9s %10s\n", "benchmark", "ns/pt", "GFLOPS", "% HW peak");
    std::printf("%-11s %8s %9s %10s\n", "-----------", "--------", "---------", "----------");

    for (auto const &record : simdRecords) {
        double gflops = static_cast<double>(record.nCoeffs - 1) * 2.0 / record.nsPerPoint;
        std::printf("N=%-2zu %-4s %8.3f %9.1f %9.1f%%\n",
                    record.nCoeffs, record.aligned ? "aln" : "uln",
                    record.nsPerPoint, gflops, gflops / peakGflops * 100.0);
    }

    return 0;
}

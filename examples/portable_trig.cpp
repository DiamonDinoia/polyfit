#include "portable_trig.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

struct TimingResult {
    std::string_view label;
    double ns_per_call = 0.0;
    double checksum = 0.0;
};

struct TimingSummary {
    std::string_view label;
    double median_ns = 0.0;
    double min_ns = 0.0;
    double max_ns = 0.0;
    double checksum = 0.0;
};

template<typename T> struct AccuracyStats {
    T max_abs_sin = T(0);
    T max_abs_cos = T(0);
    T rms_sin = T(0);
    T rms_cos = T(0);
    T worst_sin_angle = T(0);
    T worst_cos_angle = T(0);
};

template<typename T, int TolDigits = polyfit_examples::portable_trig::default_approx_digits<T>::value>
AccuracyStats<T> measure_accuracy(std::size_t samples, T range_scale) {
    AccuracyStats<T> stats;
    long double sum_sq_sin = 0.0L;
    long double sum_sq_cos = 0.0L;
    const long double pi = static_cast<long double>(polyfit_examples::portable_trig::detail::pi<T>());
    const long double range = static_cast<long double>(range_scale) * pi;

    for (std::size_t i = 0; i < samples; ++i) {
        const long double alpha = (samples > 1) ? static_cast<long double>(i) / static_cast<long double>(samples - 1) : 0.0L;
        const T angle = static_cast<T>((2.0L * alpha - 1.0L) * range);
        const std::pair<T, T> approx = polyfit_examples::portable_trig::sincos<T, TolDigits>(angle);
        const T ref_sin = std::sin(angle);
        const T ref_cos = std::cos(angle);
        const T err_sin = std::abs(approx.first - ref_sin);
        const T err_cos = std::abs(approx.second - ref_cos);

        if (err_sin > stats.max_abs_sin) {
            stats.max_abs_sin = err_sin;
            stats.worst_sin_angle = angle;
        }
        if (err_cos > stats.max_abs_cos) {
            stats.max_abs_cos = err_cos;
            stats.worst_cos_angle = angle;
        }

        sum_sq_sin += static_cast<long double>(err_sin) * static_cast<long double>(err_sin);
        sum_sq_cos += static_cast<long double>(err_cos) * static_cast<long double>(err_cos);
    }

    stats.rms_sin = static_cast<T>(std::sqrt(sum_sq_sin / static_cast<long double>(samples)));
    stats.rms_cos = static_cast<T>(std::sqrt(sum_sq_cos / static_cast<long double>(samples)));
    return stats;
}

template<typename T> std::vector<T> make_angles(std::size_t count, T range_scale) {
    std::vector<T> angles(count);
    std::uint64_t state = 0x123456789abcdefULL;
    const long double pi = static_cast<long double>(polyfit_examples::portable_trig::detail::pi<T>());
    const long double range = static_cast<long double>(range_scale) * pi;
    constexpr int unit_bits = std::numeric_limits<double>::digits;
    const long double inv_unit = 1.0L / static_cast<long double>(static_cast<std::uint64_t>(1) << unit_bits);
    constexpr int discard_bits = std::numeric_limits<std::uint64_t>::digits - unit_bits;

    for (std::size_t i = 0; i < count; ++i) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const long double unit = static_cast<long double>(state >> discard_bits) * inv_unit;
        angles[i] = static_cast<T>((2.0L * unit - 1.0L) * range);
    }
    return angles;
}

template<typename T, class Eval>
double benchmark_ns_per_call(const std::vector<T> &angles, std::size_t passes, Eval eval, T &checksum) {
    using clock = std::chrono::steady_clock;

    T sum = T(0);
    const auto start = clock::now();
    for (std::size_t pass = 0; pass < passes; ++pass) {
        for (T angle : angles) {
            const std::pair<T, T> sc = eval(angle);
            sum += sc.first * T(0.25) + sc.second * T(0.75);
        }
    }
    const auto end = clock::now();

    checksum = sum;
    const auto elapsed_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    return static_cast<double>(elapsed_ns) / static_cast<double>(angles.size() * passes);
}

template<typename T> void print_accuracy(std::string_view label, const AccuracyStats<T> &stats) {
    const std::ios::fmtflags old_flags = std::cout.flags();
    const std::streamsize old_precision = std::cout.precision();
    std::cout << std::scientific << std::setprecision(3);
    std::cout << label << " max|sin err|=" << stats.max_abs_sin << " at " << stats.worst_sin_angle
              << ", max|cos err|=" << stats.max_abs_cos << " at " << stats.worst_cos_angle
              << ", rms(sin)=" << stats.rms_sin << ", rms(cos)=" << stats.rms_cos << '\n';
    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
}

double median_of_samples(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t mid = values.size() / 2;
    if ((values.size() & 1U) != 0U) return values[mid];
    return 0.5 * (values[mid - 1] + values[mid]);
}

template<class Bench>
TimingSummary benchmark_repeated(std::string_view label, std::size_t warmup_rounds, std::size_t measured_rounds, Bench bench) {
    std::vector<double> samples;
    samples.reserve(measured_rounds);

    for (std::size_t warmup = 0; warmup < warmup_rounds; ++warmup) {
        (void)bench();
    }

    double checksum = 0.0;
    for (std::size_t round = 0; round < measured_rounds; ++round) {
        const TimingResult result = bench();
        samples.push_back(result.ns_per_call);
        checksum = result.checksum;
    }

    const auto minmax = std::minmax_element(samples.begin(), samples.end());
    return TimingSummary{label, median_of_samples(samples), *minmax.first, *minmax.second, checksum};
}

void print_timing_table(const std::vector<TimingSummary> &portable_rows, const TimingSummary &std_impl) {
    std::cout << "| implementation | median ns/call | min-max ns | checksum | vs std |\n";
    std::cout << "|---|---:|---:|---:|---:|\n";
    for (const TimingSummary &row : portable_rows) {
        const double speedup = std_impl.median_ns / row.median_ns;
        std::cout << "| " << row.label << " | " << row.median_ns << " | " << row.min_ns << "-" << row.max_ns << " | "
                  << row.checksum << " | " << speedup << "x faster |\n";
    }
    std::cout << "| std::sin + std::cos | " << std_impl.median_ns << " | " << std_impl.min_ns << "-" << std_impl.max_ns
              << " | " << std_impl.checksum << " | baseline |\n";
}

const char *feature_summary() noexcept {
#if PF_HAS_CXX23
    return "C++23 path: std::numbers pi + PF_UNREACHABLE backed by std::unreachable when available";
#elif PF_HAS_CXX20
    return "C++20 path: std::numbers pi + consteval-selected coefficient tables";
#else
    return "C++17 path: polyfit pi constant + constexpr-selected coefficient tables";
#endif
}

template<typename T>
TimingResult run_portable_benchmark(std::string_view label, const std::vector<T> &angles, std::size_t passes) {
    T checksum = T(0);
    const double ns = benchmark_ns_per_call(
        angles, passes,
        [](T angle) { return polyfit_examples::portable_trig::sincos<T>(angle); },
        checksum);
    return TimingResult{label, ns, static_cast<double>(checksum)};
}

template<typename T, int TolDigits>
TimingResult run_portable_benchmark(std::string_view label, const std::vector<T> &angles, std::size_t passes) {
    T checksum = T(0);
    const double ns = benchmark_ns_per_call(
        angles, passes,
        [](T angle) { return polyfit_examples::portable_trig::sincos<T, TolDigits>(angle); },
        checksum);
    return TimingResult{label, ns, static_cast<double>(checksum)};
}

template<typename T>
TimingResult run_std_benchmark(const std::vector<T> &angles, std::size_t passes) {
    T checksum = T(0);
    const double ns = benchmark_ns_per_call(
        angles, passes,
        [](T angle) { return std::pair<T, T>{std::sin(angle), std::cos(angle)}; },
        checksum);
    return TimingResult{"std::sin + std::cos", ns, static_cast<double>(checksum)};
}

} // namespace

int main() {
    using polyfit_examples::portable_trig::cis;
    using polyfit_examples::portable_trig::polar;
    using polyfit_examples::portable_trig::sincos;

    std::cout << std::fixed << std::setprecision(9);
    std::cout << "portable_trig example built with C++" << (PF_CPLUSPLUS / 100 % 100) << '\n';
    std::cout << feature_summary() << '\n';
    std::cout << "Build the C++17/C++20/C++23 example variants to compare standard-dependent timing.\n";
#ifndef NDEBUG
    std::cout << "Timing note: build a Release configuration for meaningful throughput numbers.\n";
#endif

    const auto unit = sincos<double>(0.5);
    const auto z1 = cis<double>(0.5);
    const auto z2 = polar<double>(2.0, 0.5);
    std::cout << "demo sincos(0.5) -> sin=" << unit.first << ", cos=" << unit.second << '\n';
    std::cout << "demo cis(0.5) -> (" << z1.real() << ", " << z1.imag() << ")\n";
    std::cout << "demo polar(2, 0.5) -> (" << z2.real() << ", " << z2.imag() << ")\n";

    const std::size_t accuracy_samples = 250001;
    std::cout << "accuracy sweeps (" << accuracy_samples << " samples each)\n";
    print_accuracy("float default   [-256*pi, 256*pi]",
                   measure_accuracy<float>(accuracy_samples, 256.0f));
    print_accuracy("double default  [-1*pi,   1*pi ]",
                   measure_accuracy<double>(accuracy_samples, 1.0));
    print_accuracy("double digits14 [-256*pi, 256*pi]",
                   measure_accuracy<double, 14>(accuracy_samples, 256.0));
    print_accuracy("double default  [-256*pi, 256*pi]",
                   measure_accuracy<double>(accuracy_samples, 256.0));
    print_accuracy("double default  [-1.000000123e6*pi, 1.000000123e6*pi]",
                   measure_accuracy<double>(accuracy_samples, 1000000.123));

    const std::size_t timing_samples = 1U << 19;
    const std::size_t timing_passes = 8;
    const auto primary_angles = make_angles<double>(timing_samples, 1.0);
    const auto moderate_angles = make_angles<double>(timing_samples, 256.0);
    constexpr std::size_t warmup_rounds = 2;
    constexpr std::size_t measured_rounds = 9;
    std::cout << "timing sweeps: " << timing_samples << " angles, " << timing_passes << " passes, "
              << warmup_rounds << " warmups, " << measured_rounds << " measured rounds\n";

    std::cout << "primary-range timing over [-pi, pi]\n";
    const std::vector<TimingSummary> primary_rows = {
        benchmark_repeated(
            "sincos<double>", warmup_rounds, measured_rounds,
            [&] {
                return run_portable_benchmark<double, 15>(
                    "sincos<double>", primary_angles, timing_passes);
            }),
        benchmark_repeated(
            "sincos<double,14>", warmup_rounds, measured_rounds,
            [&] {
                return run_portable_benchmark<double, 14>(
                    "sincos<double,14>", primary_angles, timing_passes);
            }),
    };
    const TimingSummary primary_std = benchmark_repeated(
        "std::sin + std::cos", warmup_rounds, measured_rounds,
        [&] { return run_std_benchmark(primary_angles, timing_passes); });
    print_timing_table(primary_rows, primary_std);

    std::cout << "moderate-range timing over [-256*pi, 256*pi]\n";
    const std::vector<TimingSummary> moderate_rows = {
        benchmark_repeated(
            "sincos<double,14>", warmup_rounds, measured_rounds,
            [&] {
                return run_portable_benchmark<double, 14>(
                    "sincos<double,14>", moderate_angles, timing_passes);
            }),
        benchmark_repeated(
            "sincos<double>", warmup_rounds, measured_rounds,
            [&] {
                return run_portable_benchmark<double, 15>(
                    "sincos<double>", moderate_angles, timing_passes);
            }),
    };
    const TimingSummary moderate_std = benchmark_repeated(
        "std::sin + std::cos", warmup_rounds, measured_rounds,
        [&] { return run_std_benchmark(moderate_angles, timing_passes); });
    print_timing_table(moderate_rows, moderate_std);

    return 0;
}

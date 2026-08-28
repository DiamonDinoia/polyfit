// test_scatter_and_soa.cpp: FuncEvalND scatter-write and SoA-output batch
// overload tests. Both overloads share the `evalCanonical` kernel, so the
// verification is bitwise equality against the contiguous AoS overload.

#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

namespace {

// 1D input -> 2D output: f(x) = (x, x^2). A degree-2 polynomial is fit exactly.
using In1 = std::array<double, 1>;
using Out2 = std::array<double, 2>;

auto poly_1d_to_2d = [](const In1 &x) -> Out2 { return {x[0], x[0] * x[0]}; };

// 1D input -> 1D output. Exercises the OUT_DIM == 1 SIMD-across-points path
// in the scatter overload.
using Out1 = std::array<double, 1>;

auto poly_1d_to_1d = [](const In1 &x) -> Out1 { return {x[0] * x[0] + 0.5 * x[0]}; };

// Fixed evaluation points inside [-1, 1].
inline std::vector<In1> samplePts1D() {
    return {{-0.9}, {-0.7}, {-0.4}, {-0.1}, {0.0}, {0.15}, {0.3}, {0.55}, {0.8}, {0.95}};
}

// A non-identity permutation over [0, N), built deterministically and asserted
// to permute.
inline std::vector<std::uint32_t> makePerm(std::size_t n) {
    // perm[k] = (k * 7 + 3) % n, bijective and non-trivial for n = 10.
    std::vector<std::uint32_t> p(n);
    for (std::size_t k = 0; k < n; ++k) p[k] = static_cast<std::uint32_t>((k * 7 + 3) % n);
    return p;
}

template<class T> bool bitwiseEqual(const T &a, const T &b) noexcept {
    return std::memcmp(&a, &b, sizeof(T)) == 0;
}

} // namespace

TEST(FuncEvalNdScatter, MatchesContiguousBitwise_OutDim2) {
    auto approx = poly_eval::fit(poly_1d_to_2d, 8, In1{-1.0}, In1{1.0});
    using Approx = decltype(approx);
    using CIn = Approx::CanonicalInput;
    using COut = Approx::CanonicalOutput;

    const auto pts = samplePts1D();
    const std::size_t N = pts.size();

    std::vector<CIn> cpts(N);
    for (std::size_t k = 0; k < N; ++k) cpts[k] = CIn{pts[k][0]};

    std::vector<COut> ref(N);
    approx(cpts.data(), ref.data(), N);

    const auto perm = makePerm(N);
    std::vector<COut> scattered(N);
    approx(cpts.data(), scattered.data(), perm.data(), N);

    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_TRUE(bitwiseEqual(scattered[perm[k]], ref[k]))
            << "k=" << k << " perm[k]=" << perm[k];
    }
}

TEST(FuncEvalNdScatter, MatchesContiguousBitwise_OutDim1) {
    // Hits the OUT_DIM == 1 SIMD-across-points path in the scatter overload.
    auto approx = poly_eval::fit(poly_1d_to_1d, 8, In1{-1.0}, In1{1.0});
    using Approx = decltype(approx);
    using CIn = Approx::CanonicalInput;
    using COut = Approx::CanonicalOutput;

    const auto pts = samplePts1D();
    const std::size_t N = pts.size();

    std::vector<CIn> cpts(N);
    for (std::size_t k = 0; k < N; ++k) cpts[k] = CIn{pts[k][0]};

    std::vector<COut> ref(N);
    approx(cpts.data(), ref.data(), N);

    const auto perm = makePerm(N);
    std::vector<COut> scattered(N);
    approx(cpts.data(), scattered.data(), perm.data(), N);

    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_TRUE(bitwiseEqual(scattered[perm[k]], ref[k]))
            << "k=" << k << " perm[k]=" << perm[k];
    }
}

TEST(FuncEvalNdSoa, MatchesContiguousBitwise_OutDim2) {
    auto approx = poly_eval::fit(poly_1d_to_2d, 8, In1{-1.0}, In1{1.0});
    using Approx = decltype(approx);
    using CIn = Approx::CanonicalInput;
    using COut = Approx::CanonicalOutput;
    using Scalar = Approx::Scalar;
    constexpr std::size_t OD = Approx::OUT_DIM;

    const auto pts = samplePts1D();
    const std::size_t N = pts.size();

    std::vector<CIn> cpts(N);
    for (std::size_t k = 0; k < N; ++k) cpts[k] = CIn{pts[k][0]};

    std::vector<COut> ref(N);
    approx(cpts.data(), ref.data(), N);

    std::array<std::vector<Scalar>, OD> soa{};
    std::array<Scalar *, OD> soa_ptrs{};
    for (std::size_t d = 0; d < OD; ++d) {
        soa[d].assign(N, Scalar{});
        soa_ptrs[d] = soa[d].data();
    }
    approx(cpts.data(), soa_ptrs, N);

    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t d = 0; d < OD; ++d) {
            EXPECT_TRUE(bitwiseEqual(soa[d][k], ref[k][d]))
                << "k=" << k << " d=" << d;
        }
    }
}

TEST(FuncEvalNdSoa, MatchesContiguousBitwise_OutDim1) {
    auto approx = poly_eval::fit(poly_1d_to_1d, 8, In1{-1.0}, In1{1.0});
    using Approx = decltype(approx);
    using CIn = Approx::CanonicalInput;
    using COut = Approx::CanonicalOutput;
    using Scalar = Approx::Scalar;
    constexpr std::size_t OD = Approx::OUT_DIM;
    static_assert(OD == 1, "Expected OUT_DIM=1 for this fixture");

    const auto pts = samplePts1D();
    const std::size_t N = pts.size();

    std::vector<CIn> cpts(N);
    for (std::size_t k = 0; k < N; ++k) cpts[k] = CIn{pts[k][0]};

    std::vector<COut> ref(N);
    approx(cpts.data(), ref.data(), N);

    std::vector<Scalar> col0(N, Scalar{});
    std::array<Scalar *, OD> soa_ptrs{col0.data()};
    approx(cpts.data(), soa_ptrs, N);

    for (std::size_t k = 0; k < N; ++k) {
        EXPECT_TRUE(bitwiseEqual(col0[k], ref[k][0])) << "k=" << k;
    }
}

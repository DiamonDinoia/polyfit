#include "polyfit/polyfit.hpp"
#include <gtest/gtest.h>

TEST(FuncEvalMany, PackedMatchesIndividual) {
    using namespace poly_eval;
    auto f1 = [](double x) {
        return std::sin(x);
    };
    auto f2 = [](double x) {
        return std::cos(2 * x);
    };
    double a = -1.0, b = 1.0;

    auto fe1 = fit<8>(f1, a, b);
    auto fe2 = fit<8>(f2, a, b);

    auto packed = pack(fe1, fe2);

    double x = 0.3;
    auto p = packed(x);
    auto r1 = fe1(x);
    auto r2 = fe2(x);
    EXPECT_NEAR(p[0], r1, 1e-6);
    EXPECT_NEAR(p[1], r2, 1e-6);

    // bulk evaluation with per-polynomial inputs
    std::array<double, 2> xs = {0.1, 0.2};
    auto p2 = packed(xs);
    EXPECT_NEAR(p2[0], fe1(xs[0]), 1e-6);
    EXPECT_NEAR(p2[1], fe2(xs[1]), 1e-6);
}

TEST(FuncEvalMany, EvaluateChunkMatchesBatch) {
    using namespace poly_eval;
    auto f1 = [](double x) { return std::sin(x); };
    auto f2 = [](double x) { return std::cos(2 * x); };

    auto packed = pack(fit<8>(f1, -1.0, 1.0), fit<8>(f2, -1.0, 1.0));
    using Packed = decltype(packed);

    // Sweep sizes around the chunk-block boundary plus a partial tail.
    for (std::size_t n : {std::size_t{1}, std::size_t{3}, Packed::kChunkBlockSize,
                          Packed::kChunkBlockSize * 2 + 5}) {
        std::vector<double> xs(n);
        for (std::size_t i = 0; i < n; ++i) xs[i] = -1.0 + 2.0 * double(i) / double(n + 1);

        std::vector<double> outBatch(n * Packed::COUNT);
        std::vector<double> outChunk(n * Packed::COUNT);
        packed(xs.data(), outBatch.data(), n);
        packed.evaluateChunk(xs.data(), outChunk.data(), n);
        for (std::size_t i = 0; i < outBatch.size(); ++i) EXPECT_EQ(outBatch[i], outChunk[i]);
    }

    // Drive evaluateChunk in a chunk-at-a-time loop, scatter into a permuted destination.
    const std::size_t n = Packed::kChunkBlockSize * 3 + 7;
    std::vector<double> xs(n);
    for (std::size_t i = 0; i < n; ++i) xs[i] = -1.0 + 2.0 * double(i) / double(n + 1);
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = (n - 1) - i; // simple reverse permutation

    std::vector<double> dst(n * Packed::COUNT);
    std::array<double, Packed::kChunkBlockSize * Packed::COUNT> scratch{};
    for (std::size_t off = 0; off < n; off += Packed::kChunkBlockSize) {
        const auto m = std::min(Packed::kChunkBlockSize, n - off);
        packed.evaluateChunk(xs.data() + off, scratch.data(), m);
        for (std::size_t k = 0; k < m; ++k)
            std::copy_n(scratch.begin() + static_cast<std::ptrdiff_t>(k * Packed::COUNT), Packed::COUNT,
                        dst.begin() + static_cast<std::ptrdiff_t>(perm[off + k] * Packed::COUNT));
    }

    // Reference: contiguous batch on reverse-ordered inputs.
    std::vector<double> ref(n * Packed::COUNT);
    std::vector<double> xsRev(n);
    for (std::size_t i = 0; i < n; ++i) xsRev[i] = xs[perm[i]];
    packed(xsRev.data(), ref.data(), n);
    for (std::size_t i = 0; i < dst.size(); ++i) EXPECT_EQ(dst[i], ref[i]);
}

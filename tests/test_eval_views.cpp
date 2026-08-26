// test_eval_views.cpp — external-coefficient evaluation views.
// `FuncEval` exposes its domain map (`domainInvSpan`/`domainSumEndpoints`/
// `domainIsIdentity`) and `FuncEvalND` exposes `coeffsMdspan()`/`domainView()`;
// the free `horner` / `horner_nd_batch*` entry points evaluate over those
// views. The member batch overloads forward to the same free functions, so
// the verification is bitwise equality.

#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

namespace {

template<class T> bool bitwiseEqual(const T &a, const T &b) noexcept {
    return std::memcmp(&a, &b, sizeof(T)) == 0;
}

constexpr std::size_t kNC = 8;

double cubic(double x) { return x * x * x - 2.0 * x + 1.0; }

// The bitwise-parity claim holds for the Horner scalar kernel with an unfused
// domain; the default `fit` policy (Hybrid kernel, Auto fusion) evaluates
// single points through a different kernel.
using Eval1D = poly_eval::FuncEval<double (*)(double), kNC, 1, poly_eval::FusionMode::Never,
                                   poly_eval::ScalarKernel::Horner, 0>;

} // namespace

TEST(EvalViews1D, FreeHornerBatchMatchesMemberBitwise) {
    // Non-identity domain so the map is exercised.
    constexpr double a = 0.5, b = 4.0;
    const Eval1D approx(&cubic, a, b);

    EXPECT_FALSE(approx.domainIsIdentity());

    std::vector<double> pts;
    for (int i = 0; a + 0.37 * i <= b; ++i) pts.push_back(a + 0.37 * i);
    const std::size_t n = pts.size();

    std::vector<double> ref(n), got(n);
    approx(pts.data(), ref.data(), n);

    const double invSpan = approx.domainInvSpan();
    const double sumEnds = approx.domainSumEndpoints();
    const bool   ident   = approx.domainIsIdentity();
    poly_eval::horner<kNC, false, false, 0>(pts.data(), got.data(), n, approx.coeffs().data(), kNC,
                                            [&](const auto v) {
                                                if (ident) return v;
                                                return polyfit::internal::helpers::mapFromDomainScalar(v, invSpan,
                                                                                                       sumEnds);
                                            });

    for (std::size_t i = 0; i < n; ++i) EXPECT_TRUE(bitwiseEqual(got[i], ref[i])) << "i=" << i;
}

TEST(EvalViews1D, FreeHornerSinglePointMatchesMemberBitwise) {
    constexpr double a = 0.5, b = 4.0;
    const Eval1D approx(&cubic, a, b);

    for (int i = 0; a + 0.61 * i <= b; ++i) {
        const double x   = a + 0.61 * i;
        const double ref = approx(x);
        const double xm  = approx.domainIsIdentity()
                               ? x
                               : polyfit::internal::helpers::mapFromDomainScalar(x, approx.domainInvSpan(),
                                                                                 approx.domainSumEndpoints());
        const double got = poly_eval::horner<kNC>(xm, approx.coeffs().data(), kNC);
        EXPECT_TRUE(bitwiseEqual(got, ref)) << "x=" << x;
    }
}

TEST(EvalViews1D, IdentityDomainReportsIdentity) {
    const auto approx = Eval1D([](double x) { return x * x; }, -1.0, 1.0);
    EXPECT_TRUE(approx.domainIsIdentity());
}

namespace {

using In2  = std::array<double, 2>;
using Out2 = std::array<double, 2>;

auto poly_2d_to_2d = [](const In2 &p) -> Out2 {
    return {p[0] * p[1] + 0.25 * p[0], p[0] * p[0] - p[1]};
};

std::vector<In2> samplePts2D(const In2 &lo, const In2 &hi) {
    std::vector<In2> pts;
    for (int i = 0; 0.13 * i <= 1.0; ++i)
        for (int j = 0; 0.29 * j <= 1.0; ++j) {
            const double u = 0.13 * i, v = 0.29 * j;
            pts.push_back({lo[0] + u * (hi[0] - lo[0]), lo[1] + v * (hi[1] - lo[1])});
        }
    return pts;
}

} // namespace

TEST(EvalViewsND, BatchOverViewMatchesMemberBitwise) {
    const In2 lo{0.0, -2.0}, hi{3.0, 2.0}; // non-identity domain
    const auto approx = poly_eval::fit<kNC>(poly_2d_to_2d, lo, hi);
    using Approx = decltype(approx);
    using CIn    = Approx::CanonicalInput;
    using COut   = Approx::CanonicalOutput;

    const auto pts = samplePts2D(lo, hi);
    const std::size_t n = pts.size();
    std::vector<CIn> cpts(pts.begin(), pts.end());

    std::vector<COut> ref(n), got(n);
    approx(cpts.data(), ref.data(), n);

    const auto md  = approx.coeffsMdspan();
    const auto dom = approx.domainView();
    EXPECT_FALSE(dom.identity);
    poly_eval::horner_nd_batch<kNC, 2>(md, dom, cpts.data(), got.data(), n);

    for (std::size_t i = 0; i < n; ++i) EXPECT_TRUE(bitwiseEqual(got[i], ref[i])) << "i=" << i;
}

TEST(EvalViewsND, SoaOverViewMatchesMemberBitwise) {
    const In2 lo{0.0, -2.0}, hi{3.0, 2.0};
    const auto approx = poly_eval::fit<kNC>(poly_2d_to_2d, lo, hi);
    using Approx = decltype(approx);
    using CIn    = Approx::CanonicalInput;
    using COut   = Approx::CanonicalOutput;

    const auto pts = samplePts2D(lo, hi);
    const std::size_t n = pts.size();
    std::vector<CIn> cpts(pts.begin(), pts.end());

    std::vector<COut> ref(n);
    approx(cpts.data(), ref.data(), n);

    std::array<std::vector<double>, 2> soa{std::vector<double>(n), std::vector<double>(n)};
    poly_eval::horner_nd_batch_soa<kNC, 2>(approx.coeffsMdspan(), approx.domainView(), cpts.data(),
                                           std::array<double *, 2>{soa[0].data(), soa[1].data()}, n);

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < 2; ++d)
            EXPECT_TRUE(bitwiseEqual(soa[d][i], ref[i][d])) << "i=" << i << " d=" << d;
}

TEST(EvalViewsND, ScatterOverViewMatchesMemberBitwise) {
    const In2 lo{0.0, -2.0}, hi{3.0, 2.0};
    const auto approx = poly_eval::fit<kNC>(poly_2d_to_2d, lo, hi);
    using Approx = decltype(approx);
    using CIn    = Approx::CanonicalInput;
    using COut   = Approx::CanonicalOutput;

    const auto pts = samplePts2D(lo, hi);
    const std::size_t n = pts.size();
    std::vector<CIn> cpts(pts.begin(), pts.end());

    std::vector<COut> ref(n);
    approx(cpts.data(), ref.data(), n);

    std::vector<std::uint32_t> perm(n);
    for (std::size_t k = 0; k < n; ++k) perm[k] = static_cast<std::uint32_t>((k * 7 + 3) % n);

    std::vector<COut> got(n);
    poly_eval::horner_nd_batch<kNC, 2>(approx.coeffsMdspan(), approx.domainView(), cpts.data(), got.data(),
                                       perm.data(), n);

    for (std::size_t k = 0; k < n; ++k) EXPECT_TRUE(bitwiseEqual(got[perm[k]], ref[k])) << "k=" << k;
}

TEST(EvalViewsND, RuntimeDegreeViewMatchesMemberBitwise) {
    // Runtime coefficient count: NCOEFFS=0 selects dynamic mdspan extents.
    const In2 lo{0.0, -2.0}, hi{3.0, 2.0};
    const auto approx = poly_eval::fit(poly_2d_to_2d, static_cast<int>(kNC), lo, hi);
    using Approx = decltype(approx);
    using CIn    = Approx::CanonicalInput;
    using COut   = Approx::CanonicalOutput;

    const auto pts = samplePts2D(lo, hi);
    const std::size_t n = pts.size();
    std::vector<CIn> cpts(pts.begin(), pts.end());

    std::vector<COut> ref(n), got(n);
    approx(cpts.data(), ref.data(), n);

    poly_eval::horner_nd_batch<0, 2>(approx.coeffsMdspan(), approx.domainView(), cpts.data(), got.data(), n);

    for (std::size_t i = 0; i < n; ++i) EXPECT_TRUE(bitwiseEqual(got[i], ref[i])) << "i=" << i;
}

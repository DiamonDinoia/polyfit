#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#if __cpp_lib_mdspan >= 202310L
#include <mdspan>
namespace stdex = std;
#else
#include <experimental/mdspan>
namespace stdex = std::experimental;
#endif

#include "api_types.hpp"
#include "helpers.h"
#include "macros.h"
#include "numeric_utils.h"

namespace poly_eval::detail {

namespace eft = polyfit::internal::helpers::eft;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#endif

template<std::size_t N, class X, class Y>
PF_C20CONSTEXPR Buffer<Y, N> bjorckPereyra(const Buffer<X, N> &x, const Buffer<Y, N> &y) {
    const std::size_t n = (N == 0 ? x.size() : N);
    Buffer<Y, N> a = y;

    if constexpr (std::is_arithmetic_v<Y>) {
        auto comp = makeBuffer<Y, N>(n);

        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                auto [s, es] = eft::two_sum(a[i], -a[i - 1]);
                Y subComp = es + (comp[i] - comp[i - 1]);
                Y d = Y(x[i] - x[i - k - 1]);
                Y q = s / d;
                Y r = math::fma(-q, d, s);
                comp[i] = (r + subComp) / d;
                a[i] = q;
            }
            for (std::size_t i = n - 1; i > k; --i) {
                auto [s, e] = eft::two_sum(a[i], comp[i]);
                a[i] = s;
                comp[i] = e;
            }
        }
    } else if constexpr (detail::isComplex_v<Y>) {
        using Scalar = typename Y::value_type;
        auto compRe = makeBuffer<Scalar, N>(n);
        auto compIm = makeBuffer<Scalar, N>(n);

        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                Scalar d = Scalar(x[i] - x[i - k - 1]);
                auto [sr, esr] = eft::two_sum(a[i].real(), -a[i - 1].real());
                Scalar subCompRe = esr + (compRe[i] - compRe[i - 1]);
                Scalar qr = sr / d;
                Scalar rr = math::fma(-qr, d, sr);
                compRe[i] = (rr + subCompRe) / d;

                auto [si, esi] = eft::two_sum(a[i].imag(), -a[i - 1].imag());
                Scalar subCompIm = esi + (compIm[i] - compIm[i - 1]);
                Scalar qi = si / d;
                Scalar ri = math::fma(-qi, d, si);
                compIm[i] = (ri + subCompIm) / d;

                a[i] = Y(qr, qi);
            }
            for (std::size_t i = n - 1; i > k; --i) {
                auto [sr, er] = eft::two_sum(a[i].real(), compRe[i]);
                auto [si, ei] = eft::two_sum(a[i].imag(), compIm[i]);
                a[i] = Y(sr, si);
                compRe[i] = er;
                compIm[i] = ei;
            }
        }
    } else {
        for (std::size_t k = 0; k + 1 < n; ++k) {
            for (std::size_t i = n - 1; i > k; --i) {
                a[i] = (a[i] - a[i - 1]) / static_cast<Y>(x[i] - x[i - k - 1]);
            }
        }
    }

    return a;
}

template<std::size_t N, class X, class Y>
PF_C20CONSTEXPR Buffer<Y, N> newtonToMonomial(const Buffer<Y, N> &alpha, const Buffer<X, N> &nodes) {
    const std::size_t n = alpha.size();
    constexpr std::size_t WORK_SIZE = (N == 0) ? 0 : N + 1;
    auto c = makeBuffer<Y, WORK_SIZE>(n + 1);

    if constexpr (std::is_arithmetic_v<Y>) {
        auto comp = makeBuffer<Y, WORK_SIZE>(n + 1);

        std::size_t order = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++order;
            for (std::size_t j = order; j >= 1; --j) {
                auto [p, ep] = eft::two_prod(Y(nodes[ii]), c[j]);
                Y pComp = Y(nodes[ii]) * comp[j];
                auto [s, es] = eft::two_sum(c[j - 1], -p);
                comp[j] = comp[j - 1] + es - ep - pComp;
                c[j] = s;
            }
            auto [p0, ep0] = eft::two_prod(Y(nodes[ii]), c[0]);
            Y p0Comp = Y(nodes[ii]) * comp[0];
            auto [s0, es0] = eft::two_sum(alpha[ii], -p0);
            comp[0] = es0 - ep0 - p0Comp;
            c[0] = s0;
        }
        for (std::size_t k = 0; k <= n; ++k) c[k] += comp[k];
    } else if constexpr (detail::isComplex_v<Y>) {
        using Scalar = typename Y::value_type;
        auto compRe = makeBuffer<Scalar, WORK_SIZE>(n + 1);
        auto compIm = makeBuffer<Scalar, WORK_SIZE>(n + 1);

        std::size_t order = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++order;
            Scalar node = Scalar(nodes[ii]);
            for (std::size_t j = order; j >= 1; --j) {
                auto [pr, epr] = eft::two_prod(node, c[j].real());
                Scalar prComp = node * compRe[j];
                auto [sr, esr] = eft::two_sum(c[j - 1].real(), -pr);
                compRe[j] = compRe[j - 1] + esr - epr - prComp;

                auto [pi, epi] = eft::two_prod(node, c[j].imag());
                Scalar piComp = node * compIm[j];
                auto [si, esi] = eft::two_sum(c[j - 1].imag(), -pi);
                compIm[j] = compIm[j - 1] + esi - epi - piComp;
                c[j] = Y(sr, si);
            }

            auto [pr0, epr0] = eft::two_prod(node, c[0].real());
            Scalar pr0Comp = node * compRe[0];
            auto [sr0, esr0] = eft::two_sum(alpha[ii].real(), -pr0);
            compRe[0] = esr0 - epr0 - pr0Comp;

            auto [pi0, epi0] = eft::two_prod(node, c[0].imag());
            Scalar pi0Comp = node * compIm[0];
            auto [si0, esi0] = eft::two_sum(alpha[ii].imag(), -pi0);
            compIm[0] = esi0 - epi0 - pi0Comp;

            c[0] = Y(sr0, si0);
        }
        for (std::size_t k = 0; k <= n; ++k) c[k] += Y(compRe[k], compIm[k]);
    } else {
        std::size_t order = 0;
        for (std::size_t ii = n; ii-- > 0;) {
            ++order;
            for (std::size_t j = order; j >= 1; --j) c[j] = c[j - 1] - (nodes[ii] * c[j]);
            c[0] = (-nodes[ii] * c[0]) + alpha[ii];
        }
    }

    if constexpr (N == 0) {
        c.resize(n);
        return c;
    } else {
        Buffer<Y, N> result{};
        std::copy_n(c.begin(), N, result.begin());
        return result;
    }
}

template<class X, class Y>
PF_C20CONSTEXPR std::vector<Y> bjorckPereyra(const std::vector<X> &x, const std::vector<Y> &y) {
    return bjorckPereyra<0, X, Y>(x, y);
}

template<class X, class Y>
PF_C20CONSTEXPR std::vector<Y> newtonToMonomial(const std::vector<Y> &alpha, const std::vector<X> &nodes) {
    return newtonToMonomial<0, X, Y>(alpha, nodes);
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

template<std::size_t NCOEFFS, std::size_t DIM, std::size_t OUT_DIM, std::size_t... Is>
PF_C23CONSTEVAL auto makeStaticExtentsImpl(std::index_sequence<Is...>) {
    return stdex::extents<std::size_t, ((void)Is, NCOEFFS)..., OUT_DIM>{};
}

template<std::size_t NCOEFFS, std::size_t DIM, std::size_t OUT_DIM>
using StaticExtents_t = decltype(makeStaticExtentsImpl<NCOEFFS, DIM, OUT_DIM>(std::make_index_sequence<DIM>{}));

template<class Scalar, std::size_t NCOEFFS, std::size_t DIM, std::size_t OUT_DIM>
PF_C23CONSTEVAL std::size_t storageRequired() {
    using extents_t = StaticExtents_t<NCOEFFS, DIM, OUT_DIM>;
    using mdspan_t = stdex::mdspan<Scalar, extents_t, stdex::layout_right>;
    return typename mdspan_t::mapping_type{extents_t{}}.required_span_size();
}

template<std::size_t NCOEFFS, std::size_t DIM_IN, std::size_t DIM_OUT, std::size_t... Is>
PF_C23CONSTEVAL auto makeStaticExtents(std::index_sequence<Is...>) {
    return stdex::extents<std::size_t, ((void)Is, NCOEFFS)..., DIM_OUT>{};
}

template<std::size_t COUNT = 0, typename T>
PF_C20CONSTEXPR auto linspace(const T &start, const T &end, int numPoints = COUNT) {
    Buffer<T, COUNT> points{};
    if (numPoints <= 0) {
        return points;
    }

    const auto count = static_cast<std::size_t>(numPoints);
    if constexpr (COUNT == 0) {
        points.resize(count);
    }
    if (numPoints <= 1) {
        if (numPoints == 1) points[0] = start;
        return points;
    }

    if constexpr (std::is_arithmetic_v<T>) {
        const T step = (end - start) / static_cast<T>(count - 1);
        for (std::size_t i = 0; i < count; ++i) {
            points[i] = start + (static_cast<T>(i) * step);
        }
        return points;
    } else {
        constexpr std::size_t dim = std::tuple_size_v<poly_eval::remove_cvref_t<T>>;
        using Scalar = std::remove_cv_t<decltype(start[0])>;
        T step{};
        for (std::size_t i = 0; i < dim; ++i) {
            step[i] = (end[i] - start[i]) / static_cast<Scalar>(count - 1);
        }
        for (std::size_t k = 0; k < count; ++k) {
            for (std::size_t i = 0; i < dim; ++i) {
                points[k][i] = start[i] + (static_cast<Scalar>(k) * step[i]);
            }
        }
        return points;
    }
}

template<typename T, class Step> PF_ALWAYS_INLINE constexpr void forEachComponent(const T &value, Step &&step) {
    if constexpr (detail::hasTupleSize_v<T>) {
        for (std::size_t i = 0; i < std::tuple_size_v<poly_eval::remove_cvref_t<T>>; ++i) {
            step(value[i]);
        }
    } else {
        step(value);
    }
}

template<typename T, class Step>
PF_ALWAYS_INLINE constexpr void forEachComponentPair(const T &lhs, const T &rhs, Step &&step) {
    if constexpr (detail::hasTupleSize_v<T>) {
        for (std::size_t i = 0; i < std::tuple_size_v<poly_eval::remove_cvref_t<T>>; ++i) {
            step(lhs[i], rhs[i]);
        }
    } else {
        step(lhs, rhs);
    }
}

template<typename T> PF_C20CONSTEXPR double relativeError(const T &approx, const T &actual) {
    double err = 0.0;
    forEachComponentPair(approx, actual, [&](const auto &approxValue, const auto &actualValue) {
        err = std::max(err, math::abs(1.0 - approxValue / actualValue));
    });
    return err;
}

template<typename T> PF_C20CONSTEXPR double relativeL2Norm(const T &approx, const T &actual) {
    const auto squaredNorm = [](const auto &v) constexpr noexcept -> double {
        if constexpr (detail::isComplex_v<poly_eval::remove_cvref_t<decltype(v)>>) {
            return double(v.real()) * double(v.real()) + double(v.imag()) * double(v.imag());
        } else {
            return double(v) * double(v);
        }
    };

    double numerator = 0.0;
    double denominator = 0.0;
    forEachComponentPair(approx, actual, [&](const auto &approxValue, const auto &actualValue) {
        numerator += squaredNorm(approxValue - actualValue);
        denominator += squaredNorm(actualValue);
    });

    const double ratio = denominator == 0.0 ? numerator : numerator / denominator;
    PF_IF_CONSTEVAL { return math::sqrt(ratio); }
    return std::sqrt(ratio);
}

} // namespace poly_eval::detail

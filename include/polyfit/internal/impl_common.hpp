#pragma once

// Internal: includes, namespace setup, and shared detail helpers used by
// the per-class implementation files (func_eval_impl.hpp,
// func_eval_many_impl.hpp, func_eval_nd_impl.hpp, fit_free_impl.hpp).

#include <cassert>
#include <stdexcept>
#include <xsimd/xsimd.hpp>

#include "helpers.h"
#include "macros.h"
#include "poly_eval.h"
#include "utils.h"

#include <poet/poet.hpp>

namespace poly_eval {
namespace detail {

constexpr int validatePositiveCoeffCount(const int n) {
    if (n <= 0) PF_UNLIKELY {
        throw std::invalid_argument("nCoeffs must be positive");
    }
    return n;
}

template<typename T> constexpr void validateDomain(const T &a, const T &b) {
    if constexpr (detail::isFixedIndexable_v<T>) {
        using Access = detail::FixedContainerAccess<T>;
        for (std::size_t i = 0; i < Access::size; ++i) {
            if (Access::get(a, i) == Access::get(b, i)) PF_UNLIKELY {
                throw std::invalid_argument("Domain endpoints must differ in every dimension");
            }
        }
    } else {
        if (a == b) PF_UNLIKELY {
            throw std::invalid_argument("Domain endpoints must differ");
        }
    }
}

template<class Evaluator, class Func, class Points>
[[nodiscard]] double maxRelativeError(Func F, const Evaluator &evaluator, const Points &points) {
    double maxErr = 0.0;
    for (const auto &pt : points) {
        maxErr = std::max(maxErr, relativeL2Norm(evaluator(pt), F(pt)));
    }
    return maxErr;
}

template<class Evaluator, class Func, class Spec, class InputType>
[[nodiscard]] Evaluator fitToTolerance(Func F, Spec tolerance, InputType a, InputType b, std::size_t evalPointCount,
                                       std::size_t maxCoeffCount) {
    if (tolerance <= Spec(0)) PF_UNLIKELY {
        throw std::invalid_argument("Requested error tolerance must be positive");
    }

    validateDomain(a, b);
    const auto evalPoints = linspace(a, b, int(evalPointCount));

    for (std::size_t candidateCoeffCount = 1; candidateCoeffCount <= maxCoeffCount; ++candidateCoeffCount) {
        Evaluator evaluator(F, int(candidateCoeffCount), a, b);
        if (maxRelativeError(F, evaluator, evalPoints) <= tolerance) {
            return evaluator;
        }
    }

    const Evaluator fallback(F, int(maxCoeffCount), a, b);
    throw std::runtime_error("No coefficient count found for requested error tolerance. eps=" +
                             std::to_string(tolerance) + ", maxCoeffCount=" + std::to_string(maxCoeffCount) +
                             ", maxErr=" + std::to_string(maxRelativeError(F, fallback, evalPoints)));
}

} // namespace detail
} // namespace poly_eval

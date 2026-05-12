#pragma once

#include "impl_common.hpp"
#include "func_eval_impl.hpp"
#include "func_eval_many_impl.hpp"
#include "func_eval_nd_impl.hpp"

// Out-of-class definitions for the free fit() overloads and pack().

namespace poly_eval {

template<std::size_t NCOEFFS, class Func, class... Tags>
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported or duplicate fit tag");
    static_assert(NCOEFFS > 0, "Compile-time coefficient count must be positive");
    using Evaluator = FitEvaluator<Func, NCOEFFS, Options::ITERS, Options::FUSION_MODE, ScalarKernel::Hybrid,
                                   Options::HYBRID_K>;
    return Evaluator(F, a, b);
}

template<class Func, class Spec, class... Tags>
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, Spec spec, fitInput_t<Func> a, fitInput_t<Func> b, Tags...) {
    using Options = detail::FitOptions<Tags...>;
    static_assert(Options::VALID, "Unsupported or duplicate fit tag");

    if constexpr (isIntegralLike_v<Spec>) {
        const auto nCoeffs = detail::validatePositiveCoeffCount(static_cast<int>(spec));
        if constexpr (takesNdInput_v<Func>) {
            return FuncEvalND<Func, 0, Options::FUSION_MODE, ScalarKernel::Hybrid, Options::HYBRID_K>(
                F, nCoeffs, a, b);
        } else {
            using Evaluator = FitEvaluator<Func, 0, Options::ITERS, Options::FUSION_MODE, ScalarKernel::Hybrid,
                                           Options::HYBRID_K>;
            return Evaluator(F, nCoeffs, a, b);
        }
    } else if constexpr (isFloatingPointLike_v<Spec>) {
        using Evaluator = FitEvaluator<Func, 0, Options::ITERS, Options::FUSION_MODE, ScalarKernel::Hybrid,
                                       Options::HYBRID_K>;
        return detail::fitToTolerance<Evaluator>(F, spec, a, b, Options::EVAL_POINTS, Options::MAX_NCOEFFS);
    } else {
        static_assert(alwaysFalse_v<Spec>,
                      "fit(...) expects an integral coefficient count or a floating-point tolerance");
    }
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double EPS, auto a, auto b, std::size_t MAX_NCOEFFS, std::size_t EVAL_POINTS, std::size_t ITERS, class Func>
[[nodiscard]] constexpr auto fit(Func F) {
    static_assert(MAX_NCOEFFS > 0, "Max coefficient count must be positive.");
    static_assert(EVAL_POINTS > 1, "Number of evaluation points must be greater than 1.");

    constexpr auto nCoeffs = [F] {
        constexpr auto computeError = [F](const auto &evaluator) {
            constexpr auto ep = detail::linspace<static_cast<int>(EVAL_POINTS)>(a, b);
            double maxErr = 0.0;
            for (const auto &pt : ep) {
                const auto actual = F(pt);
                const auto approx = evaluator.template operator()<false>(pt);
                maxErr = std::max(detail::relativeL2Norm(approx, actual), maxErr);
            }
            return maxErr;
        };
        int result = 0;
        poet::static_for<1, MAX_NCOEFFS + 1>([&](auto i) {
            if (result != 0) return;
            using Evaluator = std::conditional_t<takesNdInput_v<Func>, FuncEvalND<Func, i, FusionMode::Auto>,
                                                 FuncEval<Func, i, ITERS>>;
            if constexpr (computeError(Evaluator(F, a, b)) <= EPS) {
                result = i;
            }
        });
        return result;
    }();
    static_assert(nCoeffs != 0, "No coefficient count found for requested error tolerance.");
    using Evaluator = std::conditional_t<takesNdInput_v<Func>, FuncEvalND<Func, nCoeffs, FusionMode::Auto>,
                                         FuncEval<Func, nCoeffs, ITERS>>;
    return Evaluator(F, a, b);
}
#endif

template<typename... EvalTypes>
[[nodiscard]] PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...> pack(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval

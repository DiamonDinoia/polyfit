// test_ND_fusion_no_domain_params.cpp — FusionMode::Always storage layout.

#include <gtest/gtest.h>

#include <polyfit/polyfit.hpp>

#include "test_ND_fusion_shared.hpp"

using namespace polyfit_test_nd_fusion;

TEST(PolyEvalND, FusionAlwaysDoesNotStoreDomainParams) {
    using FE_fused =
        poly_eval::FuncEvalND<decltype(smooth_2d), 8, poly_eval::FusionMode::Always>;
    using FE_never =
        poly_eval::FuncEvalND<decltype(smooth_2d), 8, poly_eval::FusionMode::Never>;
    // With FusionMode::Always the per-axis domain params and identity flag
    // should collapse into the [[no_unique_address]] empty storage.
    static_assert(sizeof(FE_fused) < sizeof(FE_never),
                  "Always evaluator should be smaller than Never (no domain params stored)");
}

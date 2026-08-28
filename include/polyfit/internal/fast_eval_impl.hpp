#pragma once

// Umbrella for the per-class out-of-class implementations. Each split file is
// independently includable; this header aggregates them in dependency order.
//
// Do not include this header directly. Include `polyfit/polyfit.hpp` for the
// full fitting+evaluation API or `polyfit/polyeval.hpp` for evaluation only.

#include "impl_common.hpp"
#include "func_eval_impl.hpp"
#include "func_eval_many_impl.hpp"
#include "func_eval_nd_impl.hpp"
#include "fit_free_impl.hpp"

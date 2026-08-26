/**
 * @file polyeval.hpp
 * @brief Public header for polynomial evaluation primitives.
 *
 * This header exposes the evaluation-only API: free function templates that
 * evaluate one or many polynomials at one or many points given an array of
 * coefficients in Horner order (highest degree first). It does not include any
 * fitting machinery, so it is the right include when you already have
 * coefficients (computed offline, by another library, or by polyfit) and only
 * need to evaluate them.
 *
 * Public entry points (in namespace `poly_eval`):
 *   - `horner`           — single point or many points, 1D and ND
 *   - `hybrid`           — mixed Estrin/Horner single-point evaluator
 *   - `hybrid_transposed`,
 *     `hybrid_transposed_size`,
 *     `hybrid_transpose_coeffs` — SIMD transposed-coefficient variant
 *   - `horner_many`      — several 1D polynomials at one point
 *   - `horner_transposed` — several polynomials across many points
 *   - `horner_nd_batch`,
 *     `horner_nd_batch_soa` — across-points ND evaluation over external
 *     coefficients (`make_coeffs_mdspan` + `domain_nd_view`)
 *
 * Include `polyfit/polyfit.hpp` instead if you also need the fitting API
 * (`FuncEval`, `FuncEvalMany`, `FuncEvalND`, `fit`, `pack`).
 */

#pragma once

#include "internal/macros.hpp"
#include "internal/poly_eval.hpp"
#include "internal/horner_nd_batch.hpp"
#include "internal/macros_undef.hpp"

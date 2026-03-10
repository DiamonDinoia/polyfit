FuncEvalMany
============

``FuncEvalMany<EvalTypes...>`` packs several ``FuncEval`` instances together so their
coefficient arrays are laid out contiguously in memory, enabling SIMD-friendly batched
evaluation across all packed polynomials simultaneously.

Construction
------------

Use ``make_func_eval_many`` to pack existing evaluators:

.. code-block:: cpp

   #include <polyfit/polyfit.hpp>
   using namespace poly_eval;

   auto f1 = make_func_eval<8>([](double x){ return std::sin(x); }, 0.0, 1.0);
   auto f2 = make_func_eval<8>([](double x){ return std::cos(x); }, 0.0, 1.0);
   auto f3 = make_func_eval<8>([](double x){ return std::exp(x); }, 0.0, 1.0);

   auto packed = make_func_eval_many(f1, f2, f3);

Constraints:

- All packed evaluators must share the same sizing mode: either all compile-time coefficient count or
  all runtime coefficient count. Mixing the two in one pack is unsupported.
- All evaluators must use the same scalar type (e.g. all ``double``).

Scalar evaluation
-----------------

Evaluate all packed polynomials at the same point ``x``:

.. code-block:: cpp

   auto results = packed(0.5);
   // results[0] ≈ sin(0.5), results[1] ≈ cos(0.5), results[2] ≈ exp(0.5)

Variadic evaluation
-------------------

Evaluate each packed polynomial at its own argument:

.. code-block:: cpp

   auto results = packed(x0, x1, x2);
   // results[i] = packed_poly_i(xi)

Tuple evaluation
----------------

.. code-block:: cpp

   auto args = std::make_tuple(x0, x1, x2);
   auto results = packed(args);

Bulk SIMD evaluation
--------------------

Process many points at once using SIMD vectorization:

.. code-block:: cpp

   const std::size_t N = 1024;
   std::vector<double> xs(N), outs(N * 3);  // 3 polynomials

   packed(xs.data(), outs.data(), N);
   // outs[i*3+0] = poly0(xs[i]), outs[i*3+1] = poly1(xs[i]), outs[i*3+2] = poly2(xs[i])

The bulk overload drives SIMD Horner across all ``N`` points, amortizing loop overhead.

SIMD padding
------------

Internally, ``FuncEvalMany`` pads the coefficient arrays to the next SIMD vector width
(``kFPad``) when the number of packed polynomials is not a multiple of the SIMD width.
Padding coefficients are zero and do not affect results.

Truncation
----------

``packed.truncate(eps)`` delegates to each packed evaluator's ``truncate(eps)``, trimming
trailing near-zero coefficients from each polynomial independently:

.. code-block:: cpp

   packed.truncate(1e-14);

Runtime coefficient-count packing
---------------------------------

Runtime evaluators can also be packed. All must have the same runtime coefficient count ``nCoeffs``:

.. code-block:: cpp

   auto r1 = make_func_eval([](double x){ return std::sin(x); }, 16, 0.0, 1.0);
   auto r2 = make_func_eval([](double x){ return std::cos(x); }, 16, 0.0, 1.0);
   auto packed_rt = make_func_eval_many(r1, r2);

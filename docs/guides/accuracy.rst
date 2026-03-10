Accuracy
========

polyfit achieves near-machine-precision polynomial fitting through a combination of
compensated arithmetic in the fitting pipeline and iterative refinement during evaluation.

Compensated fitting pipeline
-----------------------------

Björck–Pereyra (Newton divided differences)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The divided-difference computation in ``bjorck_pereyra`` uses **error-free transformations
(EFT)**: every floating-point operation is accompanied by a companion correction term that
captures the rounding error exactly. These correction terms are accumulated and applied as
a compensated sum, so the final divided differences are accurate to nearly full machine
precision even at large coefficient counts (roughly 40 and above).

newton_to_monomial
^^^^^^^^^^^^^^^^^^

Converting Newton-form coefficients to the monomial basis involves a triangular recurrence.
``newton_to_monomial`` applies the same EFT + compensated summation strategy to prevent
error accumulation across the recurrence steps.

Chebyshev node stability
^^^^^^^^^^^^^^^^^^^^^^^^^

Chebyshev nodes are computed via a Cody–Waite minimax approximation of cosine
(``detail::cos``) rather than ``std::cos``. This ensures bit-identical node positions
across C++17 and C++20 (and across MSVC vs GCC/Clang), which matters for large fits
where divided differences become sensitive to small perturbations in node placement.

.. note::

   Using ``std::cos`` for node placement can cause near-degenerate divided differences at
   coefficient counts around 44 and above, because ``std::cos`` and the Cody–Waite approximation differ by up to 1 ULP.
   The library always uses the Cody–Waite path to avoid this instability.

Iterative refinement
--------------------

After the initial fit, ``FuncEval`` applies ``RefineIters`` rounds of iterative
refinement. Each round:

1. Evaluates the current polynomial at the Chebyshev nodes using **compensated Horner**
   (Kahan-style accumulation to avoid catastrophic cancellation in the evaluation step).
2. Computes the residual ``f(x_i) - p(x_i)`` at each node.
3. Fits a correction polynomial to the residuals.
4. Adds the correction to the coefficient array.

Compensated Horner is critical in step 1: plain Horner at large coefficient counts can itself introduce
rounding errors larger than the residual being corrected, causing the refinement to diverge.

Recommended refinement counts:

+----------+---------------+
| Coeffs   | ``Iters``     |
+==========+===============+
| ≤ 40     | 1 (default)   |
+----------+---------------+
| 40 – 80  | 2             |
+----------+---------------+
| 80 – 128 | 3             |
+----------+---------------+

Truncation
----------

``FuncEval::truncate(eps)`` removes leading coefficients whose absolute value is ≤ ``eps``,
reducing the effective coefficient count and evaluation cost without measurably affecting
accuracy (since those coefficients contribute less than ``eps`` to the output).

The adaptive overload of ``make_func_eval`` (epsilon-based) returns the smallest coefficient count that
passes the requested error check; it does not apply automatic truncation afterward.

C++ standard and platform portability
--------------------------------------

- **C++17**: Fully supported runtime fitting and evaluation. Cody–Waite cosine is always
  used for Chebyshev nodes.
- **C++20**: Adds constexpr fitting paths. Node computation uses the same Cody–Waite
  approximation as C++17, ensuring identical results across standards.
- **MSVC**: Requires ``/Zc:__cplusplus`` for C++20 constexpr paths. All fitting arithmetic
  is compatible with MSVC's floating-point model.

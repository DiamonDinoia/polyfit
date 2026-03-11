pack
====

``poly_eval::pack(...)`` combines several 1D ``FuncEval`` objects into one ``FuncEvalMany``.

Basic use
---------

.. code-block:: cpp

   auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
   auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);

   auto packed = poly_eval::pack(sinApprox, cosApprox);
   auto y = packed(0.5);

Supported call forms
--------------------

- ``packed(x)`` evaluates every polynomial at the same point
- ``packed(xs)`` evaluates each polynomial at its corresponding point
- ``packed(pts, out, count)`` evaluates many points into a flat output buffer

Constraints
-----------

- ``pack(...)`` is for 1D evaluators.
- All packed evaluators must use the same sizing mode.
- Runtime-sized evaluators may be packed together, but they must agree on the active coefficient count.

Accuracy
========

polyfit samples on Chebyshev nodes, solves in Newton form, converts to monomial coefficients, and evaluates with Horner's method.

Practical notes
---------------

- Higher-degree fits may benefit from ``poly_eval::Iters<N>``.
- The adaptive overload checks the fit on a finite grid controlled by ``EvalPts<N>``.
- ``FuseAuto`` is the default domain-mapping mode for 1D fits.

Truncation
----------

``FuncEval::truncate(eps)`` removes leading small coefficients on runtime-sized 1D fits.

``FuncEvalMany::truncate(eps)`` trims trailing coefficient rows across the packed evaluators.

These operations reduce work, but they are optional and explicit.

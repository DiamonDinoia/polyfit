Installation
============

Requirements
------------

- C++17-conforming compiler (GCC ≥ 10, Clang ≥ 12, MSVC 2019+, AppleClang ≥ 13)
- CMake ≥ 3.14
- (Optional) C++20 for constexpr fitting paths
- (Optional) C++23/C++26 for newer compile-time and constexpr-`<cmath>` paths

Feature levels
--------------

Use this as the short compatibility guide:

.. list-table::
   :header-rows: 1

   * - Toolchain level
     - What it enables
   * - C++17
     - Runtime 1D fitting, runtime ND fitting, packed 1D evaluation
   * - C++20
     - Constant-evaluated fixed-count 1D and fixed-count ND fits
   * - ``PF_HAS_CONSTEXPR_EPS_OVERLOAD``
     - Compile-time epsilon-driven 1D fitting
   * - C++23 / newer stdlib support
     - Native ``if consteval`` and cleaner compile-time dispatch
   * - C++26-capable stdlib
     - Native constexpr ``<cmath>`` paths

Tested support
--------------

The current CI matrix exercises:

.. list-table::
   :header-rows: 1

   * - Area
     - Current CI coverage
   * - Runtime 1D, runtime ND, ``pack(...)``
     - Linux ``gcc``, ``gcc-13``, ``gcc-14``, ``llvm``, ``llvm-18``, ``llvm-21``; macOS ``apple-clang``; Windows ``msvc``
   * - Fixed-count ``constexpr`` 1D and ND
     - Built and tested in the Linux C++20 and C++23 matrix
   * - Compile-time epsilon 1D
     - Exercised only when ``PF_HAS_CONSTEXPR_EPS_OVERLOAD`` is true; in current CI that means the Linux GCC jobs

polyfit is header-only. The following methods all make ``polyfit::polyfit`` available as a CMake
INTERFACE target.

FetchContent
------------

.. code-block:: cmake

   include(FetchContent)
   FetchContent_Declare(polyfit
     GIT_REPOSITORY https://github.com/DiamonDinoia/polyfit.git
     GIT_TAG        main
   )
   FetchContent_MakeAvailable(polyfit)
   target_link_libraries(my_target PRIVATE polyfit::polyfit)

CPM.cmake
---------

.. code-block:: cmake

   CPMAddPackage(
     NAME polyfit
     GITHUB_REPOSITORY DiamonDinoia/polyfit
     GIT_TAG main
   )
   target_link_libraries(my_target PRIVATE polyfit::polyfit)

find_package (after install)
----------------------------

.. code-block:: bash

   # Install once:
   cmake -S . -B build && cmake --install build --prefix /usr/local

.. code-block:: cmake

   # In consuming project:
   find_package(polyfit REQUIRED)
   target_link_libraries(my_target PRIVATE polyfit::polyfit)

Public headers
--------------

The supported public entry points are:

* ``include/polyfit/polyfit.hpp``: umbrella include for the full public API
* ``include/polyfit/polyeval.hpp``: direct include for the polynomial-evaluation API

Internal headers remain under ``include/polyfit/internal`` and are not intended for direct use.

Generated amalgamation
----------------------

If you want a generated single-header artifact for vendoring or inspection, generate it with:

.. code-block:: bash

   python3 scripts/amalgamate.py   # writes include/polyfit/internal/polyfit_amalgamated.hpp

Or build the CMake target:

.. code-block:: bash

   cmake --build build --target amalgamate

For normal use, include one of the public headers:

.. code-block:: cpp

   #include <polyfit/polyfit.hpp>

.. warning::

   The single header **does not bundle** its external dependencies. You must have the following
   libraries on your compiler include path:

   * `xsimd <https://github.com/xtensor-stack/xsimd>`_ ≥ 13.0 — SIMD abstraction layer
   * `poet <https://github.com/DiamonDinoia/poet>`_ — parallel iteration utilities
   * `kokkos/mdspan <https://github.com/kokkos/mdspan>`_ ≥ 0.6 — ``std::mdspan`` (C++17/20)
     *or* a compiler with native ``std::mdspan`` (C++23+)

   Example compile command:

   .. code-block:: bash

      g++ -std=c++17 -I path/to/xsimd/include -I path/to/poet/include \
          -I path/to/mdspan/include -I path/to/polyfit my_program.cpp

Header copy
-----------

Copy the full ``include/polyfit/`` directory into your project. You will also need
`xsimd <https://github.com/xtensor-stack/xsimd>`_, `poet <https://github.com/DiamonDinoia/poet>`_,
and (for C++17) `kokkos/mdspan <https://github.com/kokkos/mdspan>`_ on your include path.

MSVC notes
----------

For C++20 constexpr fitting on MSVC, add ``/Zc:__cplusplus`` to your compile flags so that
MSVC correctly reports the language standard version:

.. code-block:: cmake

   target_compile_options(my_target PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>)

Optional modern paths
---------------------

polyfit keeps the documented C++17 baseline and C++20 constexpr-fit support.
When newer toolchains are available, the headers enable additional paths
automatically:

* C++23: native ``if consteval`` and constexpr local statics where supported
* C++26-capable standard libraries: native constexpr ``<cmath>`` functions

These are compatibility improvements only. You do not need to opt in through a
separate CMake option.

Building the docs
-----------------

Build the docs through the CMake ``docs`` target:

.. code-block:: bash

   cmake -S . -B build-docs-check
   python3 -m pip install -r docs/requirements.txt
   cmake --build build-docs-check --target docs

This target runs the configured Doxygen and Sphinx steps for the build tree.
The generated HTML is written under ``build-docs-check/docs/html/``.

Running benchmarks
------------------

Build and run the benchmark suite through the CMake benchmark target:

.. code-block:: bash

   cmake -S . -B build-bench -DPOLYFIT_BUILD_TESTS=ON -DPOLYFIT_BUILD_EXAMPLES=OFF -DPOLYFIT_ENABLE_BENCH_NATIVE_TUNING=ON
   cmake --build build-bench --target run_benchmarks

For stable measurements, use a tuned machine state: prefer a ``performance``
CPU governor, run ``pyperf system tune`` when available, and avoid
sanitizer-enabled builds.

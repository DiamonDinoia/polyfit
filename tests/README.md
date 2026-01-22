Tests and SIMD variants
=======================

This folder contains the test CMake entry points and a small helper script to run tests locally across multiple CPU ISA variants.

```bash
cmake -S . -B build -DMONOFIT_BUILD_TESTS=ON -DMONOFIT_BUILD_SIMD_VARIANTS=ON
cmake --build build -j
ctest --test-dir build
```

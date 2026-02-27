Tests and SIMD variants
=======================

This folder contains the test CMake entry points and a small helper script to run tests locally across multiple CPU ISA variants.

```bash
cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build
```

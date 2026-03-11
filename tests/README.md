Tests
=====

Build the test targets before running ``ctest``.

```bash
cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

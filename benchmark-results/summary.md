# Benchmark Summary

All values are **ns/op** (lower is better).


| Benchmark | gcc-14 | gcc-15 | llvm-20 | llvm-21 |
| --- | ---: | ---: | ---: | ---: |
| `bench1D_cxx17::float constexpr fit` | 447.69 | 436.62 | 745.94 | 722.11 |
| `bench1D_cxx17::float coefficient-count fit` | 1001.72 | 948.14 | 956.64 | 947.49 |
| `bench1D_cxx17::float eps fit` | 1838.71 | 2301.00 | 2438.80 | 2229.96 |
| `bench1D_cxx17::float eval` | 2.11 | 1.76 | 1.76 | 2.09 |
| `bench1D_cxx17::float eval_many` | 0.29 | 0.26 | 0.17 | 0.17 |
| `bench1D_cxx17::double constexpr fit` | 468.06 | 486.68 | 791.96 | 749.10 |
| `bench1D_cxx17::double coefficient-count fit` | 1146.18 | 1089.51 | 1048.07 | 1031.85 |
| `bench1D_cxx17::double eps fit` | 1988.74 | 2414.37 | 2393.62 | 2261.21 |
| `bench1D_cxx17::double eval` | 2.11 | 2.05 | 1.78 | 2.08 |
| `bench1D_cxx17::double eval_many` | 0.55 | 0.50 | 0.34 | 0.37 |
| `bench1D_cxx17::complex<double> constexpr fit` | 610.43 | 585.42 | 674.78 | 630.26 |
| `bench1D_cxx17::complex<double> coefficient-count fit` | 1372.69 | 1389.68 | 1534.65 | 1470.27 |
| `bench1D_cxx17::complex<double> eps fit` | 2212.08 | 2854.60 | 2709.82 | 2574.88 |
| `bench1D_cxx17::complex<double> eval` | 2.31 | 2.29 | 1.76 | 2.09 |
| `bench1D_cxx17::complex<double> eval_many` | 0.87 | 1.35 | 1.19 | 1.38 |
| `bench1D_cxx17::complex<float> constexpr fit` | 632.20 | 580.79 | 926.85 | 851.36 |
| `bench1D_cxx17::complex<float> coefficient-count fit` | 1355.98 | 1270.68 | 1664.90 | 1576.94 |
| `bench1D_cxx17::complex<float> eps fit` | 2301.83 | 2670.29 | 2800.73 | 2695.56 |
| `bench1D_cxx17::complex<float> eval` | 3.60 | 3.87 | 3.47 | 3.49 |
| `bench1D_cxx17::complex<float> eval_many` | 0.49 | 0.69 | 0.59 | 0.65 |
| `bench1D_cxx20::float constexpr fit` | 438.16 | 437.07 | 748.93 | 713.79 |
| `bench1D_cxx20::float coefficient-count fit` | 995.00 | 998.16 | 923.48 | 942.59 |
| `bench1D_cxx20::float eps fit` | 1834.03 | 2331.83 | 2373.68 | 2254.94 |
| `bench1D_cxx20::float eval` | 2.11 | 1.81 | 1.82 | 2.08 |
| `bench1D_cxx20::float eval_many` | 0.30 | 0.26 | 0.17 | 0.17 |
| `bench1D_cxx20::double constexpr fit` | 454.59 | 476.26 | 810.25 | 741.60 |
| `bench1D_cxx20::double coefficient-count fit` | 1139.08 | 1175.73 | 1042.86 | 1035.06 |
| `bench1D_cxx20::double eps fit` | 1999.43 | 2505.84 | 2397.46 | 2255.26 |
| `bench1D_cxx20::double eval` | 2.11 | 2.02 | 1.79 | 2.09 |
| `bench1D_cxx20::double eval_many` | 0.55 | 0.50 | 0.34 | 0.34 |
| `bench1D_cxx20::complex<double> constexpr fit` | 582.66 | 582.71 | 674.16 | 637.13 |
| `bench1D_cxx20::complex<double> coefficient-count fit` | 1285.89 | 1380.00 | 1530.82 | 1450.21 |
| `bench1D_cxx20::complex<double> eps fit` | 2213.76 | 2936.85 | 2703.85 | 2574.77 |
| `bench1D_cxx20::complex<double> eval` | 2.26 | 2.33 | 1.77 | 2.07 |
| `bench1D_cxx20::complex<double> eval_many` | 0.88 | 1.35 | 1.32 | 1.39 |
| `bench1D_cxx20::complex<float> constexpr fit` | 621.42 | 580.01 | 927.69 | 849.09 |
| `bench1D_cxx20::complex<float> coefficient-count fit` | 1394.54 | 1333.18 | 1647.94 | 1562.89 |
| `bench1D_cxx20::complex<float> eps fit` | 2352.64 | 2772.02 | 2740.65 | 2640.44 |
| `bench1D_cxx20::complex<float> eval` | 3.59 | 3.86 | 3.51 | 3.55 |
| `bench1D_cxx20::complex<float> eval_many` | 0.49 | 0.68 | 0.67 | 0.66 |
| `benchMany_cxx17::1 funcs` | 7.22 | 7.64 | 5.31 | 5.65 |
| `benchMany_cxx17::1 funcs (non-many)` | 4.82 | 4.40 | 4.58 | 5.15 |
| `benchMany_cxx17::2 funcs` | 8.68 | 7.94 | 5.33 | 5.87 |
| `benchMany_cxx17::2 funcs (non-many)` | 9.63 | 9.30 | 8.56 | 10.00 |
| `benchMany_cxx17::3 funcs` | 8.28 | 8.61 | 5.12 | 6.14 |
| `benchMany_cxx17::3 funcs (non-many)` | 19.88 | 13.43 | 12.87 | 15.09 |
| `benchMany_cxx17::4 funcs` | 8.18 | 8.44 | 5.33 | 5.65 |
| `benchMany_cxx17::4 funcs (non-many)` | 18.95 | 17.96 | 19.38 | 19.36 |
| `benchMany_cxx17::5 funcs` | 9.01 | 10.93 | 7.92 | 8.13 |
| `benchMany_cxx17::5 funcs (non-many)` | 23.64 | 21.56 | 21.30 | 24.07 |
| `benchMany_cxx17::6 funcs` | 9.17 | 10.99 | 8.12 | 8.09 |
| `benchMany_cxx17::6 funcs (non-many)` | 28.54 | 38.16 | 25.61 | 29.44 |
| `benchMany_cxx17::7 funcs` | 9.59 | 11.80 | 8.30 | 8.42 |
| `benchMany_cxx17::7 funcs (non-many)` | 31.27 | 32.81 | 29.72 | 33.41 |
| `benchMany_cxx17::8 funcs` | 8.94 | 11.02 | 8.15 | 8.13 |
| `benchMany_cxx17::8 funcs (non-many)` | 37.46 | 34.93 | 33.74 | 38.43 |
| `benchMany_cxx17::9 funcs` | 12.21 | 14.50 | 12.68 | 9.57 |
| `benchMany_cxx17::9 funcs (non-many)` | 59.80 | 39.08 | 38.10 | 43.28 |
| `benchMany_cxx17::10 funcs` | 11.84 | 14.86 | 13.10 | 11.80 |
| `benchMany_cxx17::10 funcs (non-many)` | 52.14 | 43.62 | 42.70 | 47.94 |
| `benchMany_cxx17::11 funcs` | 14.20 | 16.22 | 11.07 | 11.21 |
| `benchMany_cxx17::11 funcs (non-many)` | 48.83 | 62.41 | 46.61 | 52.35 |
| `benchMany_cxx17::12 funcs` | 13.08 | 15.14 | 12.75 | 11.72 |
| `benchMany_cxx17::12 funcs (non-many)` | 79.70 | 52.44 | 52.63 | 61.92 |
| `benchMany_cxx17::13 funcs` | 15.93 | 16.85 | 14.21 | 13.81 |
| `benchMany_cxx17::13 funcs (non-many)` | 81.42 | 73.88 | 56.97 | 67.07 |
| `benchMany_cxx17::14 funcs` | 15.21 | 17.51 | 16.66 | 14.46 |
| `benchMany_cxx17::14 funcs (non-many)` | 69.41 | 65.19 | 61.40 | 72.13 |
| `benchMany_cxx17::15 funcs` | 15.96 | 18.90 | 15.35 | 13.61 |
| `benchMany_cxx17::15 funcs (non-many)` | 66.39 | 88.65 | 67.76 | 82.03 |
| `benchMany_cxx17::16 funcs` | 18.62 | 16.85 | 16.56 | 13.52 |
| `benchMany_cxx17::16 funcs (non-many)` | 100.43 | 73.13 | 69.86 | 88.27 |
| `benchMany_cxx17::8 funcs (mixed domains)` | 8.94 | 11.04 | 7.92 | 8.39 |
| `benchMany_cxx17::8 funcs (mixed domains, non-many)` | 37.87 | 34.99 | 33.78 | 38.52 |
| `benchMany_cxx17::16 funcs (mixed domains)` | 15.58 | 16.88 | 17.07 | 13.80 |
| `benchMany_cxx17::16 funcs (mixed domains, non-many)` | 100.41 | 73.14 | 69.91 | 89.05 |
| `benchMany_cxx17::8 funcs (complex outputs, bulk)` | 26.04 | 38.60 | 39.57 | 43.57 |
| `benchMany_cxx20::1 funcs` | 6.61 | 7.53 | 5.28 | 5.65 |
| `benchMany_cxx20::1 funcs (non-many)` | 4.82 | 4.37 | 4.57 | 5.14 |
| `benchMany_cxx20::2 funcs` | 7.81 | 7.86 | 5.29 | 5.84 |
| `benchMany_cxx20::2 funcs (non-many)` | 9.65 | 9.28 | 8.61 | 9.98 |
| `benchMany_cxx20::3 funcs` | 7.37 | 8.49 | 5.06 | 6.14 |
| `benchMany_cxx20::3 funcs (non-many)` | 19.88 | 13.46 | 12.87 | 15.03 |
| `benchMany_cxx20::4 funcs` | 7.75 | 8.81 | 5.28 | 5.65 |
| `benchMany_cxx20::4 funcs (non-many)` | 19.15 | 18.04 | 16.96 | 19.37 |
| `benchMany_cxx20::5 funcs` | 9.32 | 10.95 | 8.20 | 8.22 |
| `benchMany_cxx20::5 funcs (non-many)` | 23.62 | 21.63 | 21.33 | 24.03 |
| `benchMany_cxx20::6 funcs` | 9.36 | 10.97 | 8.47 | 8.32 |
| `benchMany_cxx20::6 funcs (non-many)` | 29.08 | 38.22 | 25.60 | 29.40 |
| `benchMany_cxx20::7 funcs` | 10.03 | 11.79 | 8.21 | 8.37 |
| `benchMany_cxx20::7 funcs (non-many)` | 31.30 | 32.84 | 29.74 | 35.00 |
| `benchMany_cxx20::8 funcs` | 9.29 | 11.05 | 8.04 | 8.36 |
| `benchMany_cxx20::8 funcs (non-many)` | 37.75 | 35.01 | 33.74 | 38.48 |
| `benchMany_cxx20::9 funcs` | 12.24 | 14.31 | 13.05 | 11.10 |
| `benchMany_cxx20::9 funcs (non-many)` | 59.81 | 39.05 | 38.12 | 42.93 |
| `benchMany_cxx20::10 funcs` | 11.84 | 14.67 | 13.59 | 9.94 |
| `benchMany_cxx20::10 funcs (non-many)` | 52.03 | 43.61 | 42.14 | 47.82 |
| `benchMany_cxx20::11 funcs` | 14.68 | 16.24 | 11.57 | 11.51 |
| `benchMany_cxx20::11 funcs (non-many)` | 48.83 | 62.45 | 46.59 | 52.63 |
| `benchMany_cxx20::12 funcs` | 13.45 | 15.14 | 12.82 | 11.25 |
| `benchMany_cxx20::12 funcs (non-many)` | 79.75 | 52.48 | 54.54 | 62.79 |
| `benchMany_cxx20::13 funcs` | 16.46 | 16.90 | 14.22 | 13.57 |
| `benchMany_cxx20:::wavy_dash: `13 funcs (non-many)` (Unstable with ~277.3 iters. Increase `minEpochIterations` to e.g. 2773)` | 81.53 | — | — | — |
| `benchMany_cxx20::14 funcs` | 15.20 | 17.31 | 16.72 | 13.76 |
| `benchMany_cxx20::14 funcs (non-many)` | 69.48 | 65.20 | 63.62 | 73.20 |
| `benchMany_cxx20::15 funcs` | 17.10 | 18.93 | 14.82 | 12.74 |
| `benchMany_cxx20::15 funcs (non-many)` | 66.28 | 88.61 | 65.81 | 82.90 |
| `benchMany_cxx20::16 funcs` | 17.88 | 16.89 | 16.40 | 33.92 |
| `benchMany_cxx20::16 funcs (non-many)` | 100.55 | 73.10 | 69.99 | 88.18 |
| `benchMany_cxx20::8 funcs (mixed domains)` | 8.91 | 11.04 | 8.39 | 8.15 |
| `benchMany_cxx20::8 funcs (mixed domains, non-many)` | 37.80 | 34.97 | 33.82 | 38.42 |
| `benchMany_cxx20::16 funcs (mixed domains)` | 15.53 | 16.94 | 15.75 | 13.89 |
| `benchMany_cxx20::16 funcs (mixed domains, non-many)` | 100.56 | 73.26 | 69.84 | 88.18 |
| `benchMany_cxx20::8 funcs (complex outputs, bulk)` | 28.34 | 38.69 | 39.45 | 43.29 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 build` | 53459.58 | 60594.58 | 82121.82 | 78378.34 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16 build` | 79275.64 | 89153.42 | 106860.79 | 113133.87 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16 build` | 1355227.25 | 1508060.29 | 2024506.80 | 1954344.50 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8 build` | 212781.43 | 235895.69 | 388682.67 | 345912.17 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8 build` | 280087.78 | 312135.41 | 521013.36 | 474820.08 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8 build` | 2221914.60 | 2555243.20 | 4234456.00 | 3850991.67 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8 build` | 2934429.75 | 3372605.00 | 5903776.50 | 5137887.50 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16` | 48.43 | 52.26 | 1301.79 | 902.23 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16` | 94.99 | 112.36 | 1240.59 | 1005.50 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16` | 1137.76 | 1289.01 | 17986.24 | 14808.66 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8` | 178.71 | 197.15 | 2824.77 | 2135.23 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8` | 94.13 | 104.67 | 2161.25 | 1962.52 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8` | 2465.87 | 2603.93 | 21473.68 | 17228.07 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8` | 1296.13 | 1719.58 | 19336.04 | 16162.50 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic point` | 48.35 | 52.27 | 1304.06 | 900.81 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical point` | 48.35 | 52.26 | 1301.93 | 901.27 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 array` | 48.38 | 52.27 | 1297.80 | 900.74 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 variadic` | 78.46 | 92.87 | 1299.10 | 900.52 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 manual loop` | 49025.16 | 52840.90 | 1339004.90 | 923429.67 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical batch` | 95093.55 | 111623.91 | 1304914.88 | 924183.17 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic loop` | 49083.57 | 52838.93 | 1339381.40 | 921979.23 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 container batch` | 93781.19 | 113605.43 | 1333358.27 | 922067.12 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 build` | 53443.01 | 61050.76 | 82517.30 | 78592.85 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16 build` | 79303.74 | 88770.29 | 107172.23 | 106271.33 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16 build` | 1348696.22 | 1518223.87 | 2082905.83 | 1967188.33 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8 build` | 231046.89 | 232864.91 | 390301.03 | 348725.73 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8 build` | 303119.50 | 307947.35 | 514868.74 | 470995.12 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8 build` | 2266123.40 | 2540415.20 | 4438255.50 | 3870364.33 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8 build` | 2995210.67 | 3348587.67 | 5560289.50 | 5042165.00 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16` | 48.43 | 52.24 | 1304.29 | 901.33 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16` | 95.07 | 112.19 | 1252.51 | 997.96 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16` | 1138.14 | 1266.69 | 18415.31 | 14821.25 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8` | 178.72 | 197.32 | 2763.07 | 2184.10 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8` | 94.09 | 107.30 | 2191.76 | 1979.52 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8` | 2468.96 | 2601.35 | 20931.97 | 17254.42 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8` | 1294.03 | 1718.67 | 19425.35 | 16133.41 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic point` | 48.44 | 52.23 | 1298.47 | 900.52 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical point` | 48.40 | 52.26 | 1299.43 | 900.19 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 array` | 48.42 | 52.22 | 1304.03 | 900.10 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 variadic` | 76.97 | 97.22 | 1303.72 | 901.56 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 manual loop` | 49011.13 | 52905.15 | 1336029.06 | 921958.48 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical batch` | 93170.48 | 112629.63 | 1305522.06 | 924962.09 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 span batch` | 76550.91 | 92145.19 | 1307997.33 | 922053.91 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic loop` | 49040.92 | 52830.85 | 1336859.76 | 922899.84 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 container batch` | 94487.65 | 114679.34 | 1330551.75 | 920658.09 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.40 | 2.18 | 3.83 | 3.42 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.28 | 0.35 | 0.34 | 0.34 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.27 | 0.36 | 0.34 | 0.35 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.55 | 4.39 | 4.62 | 5.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.61 | 0.67 | 0.71 | 0.68 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.59 | 0.67 | 0.71 | 0.68 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.55 | 7.55 | 6.73 | 7.66 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 0.95 | 1.05 | 1.13 | 1.06 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 0.93 | 1.05 | 1.14 | 1.06 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 10.18 | 11.18 | 10.21 | 10.57 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.28 | 1.45 | 1.56 | 1.43 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.26 | 1.44 | 1.56 | 1.42 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.66 | 2.45 | 5.59 | 5.31 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.84 | 2.53 | 4.25 | 4.60 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.67 | 2.29 | 3.88 | 4.03 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.33 | 2.10 | 3.65 | 3.85 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 4.88 | 4.65 | 7.31 | 8.35 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 4.76 | 4.66 | 5.79 | 7.61 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.60 | 4.50 | 5.16 | 6.46 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.32 | 4.14 | 4.85 | 6.08 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 6.75 | 8.00 | 10.31 | 11.46 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 7.59 | 7.67 | 8.39 | 10.60 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 7.01 | 7.49 | 7.67 | 9.56 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 6.92 | 7.11 | 7.43 | 9.11 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 9.76 | 12.72 | 12.74 | 14.80 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 10.74 | 11.87 | 11.80 | 13.67 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 10.19 | 11.51 | 10.72 | 12.64 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 9.98 | 10.58 | 10.92 | 12.24 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 7.42 | 8.79 | 6.39 | 6.09 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.26 | 3.07 | 2.92 | 2.73 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 2.12 | 2.36 | 3.20 | 2.90 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 1.87 | 2.03 | 2.21 | 2.02 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 14.33 | 18.67 | 13.98 | 12.86 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 4.87 | 6.65 | 7.17 | 6.44 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 18.73 | 4.25 | 6.01 | 5.36 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 3.75 | 4.32 | 4.58 | 4.47 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 21.24 | 27.99 | 21.06 | 20.30 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 7.75 | 9.91 | 11.05 | 8.94 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 5.86 | 8.27 | 10.03 | 9.28 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 5.67 | 6.74 | 7.44 | 6.90 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 28.15 | 37.33 | 30.21 | 29.66 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 10.74 | 14.32 | 13.30 | 13.20 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 8.30 | 11.83 | 13.79 | 12.39 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 7.54 | 9.08 | 10.46 | 9.58 |
| `bench_horner_cxx17::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 8.89 | 8.09 | 7.23 | 7.48 |
| `bench_horner_cxx17::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 36.14 | 27.37 | 16.90 | 16.59 |
| `bench_horner_cxx17::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 119.68 | 132.30 | 61.40 | 63.52 |
| `bench_horner_cxx17::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 265.35 | 332.31 | 167.78 | 165.33 |
| `bench_horner_cxx17::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 38.13 | 35.84 | 34.33 | 36.09 |
| `bench_horner_cxx17::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 296.99 | 233.20 | 205.77 | 225.81 |
| `bench_horner_cxx17::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2440.62 | 2127.52 | 1797.92 | 1932.82 |
| `bench_horner_cxx17::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 151.85 | 137.57 | 128.30 | 129.05 |
| `bench_horner_cxx17::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 2355.06 | 2069.84 | 1563.99 | 1593.22 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.39 | 2.17 | 3.52 | 3.42 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.28 | 0.35 | 0.34 | 0.34 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.27 | 0.36 | 0.34 | 0.35 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.55 | 4.37 | 4.66 | 5.28 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.62 | 0.67 | 0.71 | 0.68 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.59 | 0.67 | 0.71 | 0.68 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.54 | 7.58 | 6.74 | 7.72 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 0.97 | 1.05 | 1.14 | 1.05 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 0.93 | 1.05 | 1.13 | 1.05 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 10.18 | 11.18 | 10.20 | 10.84 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.31 | 1.44 | 1.56 | 1.43 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.26 | 1.43 | 1.56 | 1.43 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.65 | 2.63 | 5.71 | 5.32 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.92 | 2.55 | 4.30 | 4.61 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.69 | 2.31 | 3.85 | 4.03 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.42 | 2.06 | 3.65 | 3.84 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 5.02 | 4.66 | 8.15 | 8.33 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 5.10 | 4.67 | 5.78 | 7.55 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.57 | 4.50 | 5.21 | 6.48 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.33 | 4.23 | 4.87 | 6.08 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 6.72 | 8.00 | 10.60 | 11.58 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 7.56 | 7.68 | 8.37 | 10.54 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 7.01 | 7.50 | 7.75 | 9.53 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 6.85 | 7.12 | 7.42 | 9.10 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 9.85 | 12.72 | 13.63 | 14.81 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 10.67 | 11.87 | 11.92 | 13.66 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 10.28 | 11.51 | 10.72 | 12.62 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 9.99 | 10.58 | 10.94 | 12.24 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 7.43 | 8.79 | 6.18 | 6.09 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.26 | 3.08 | 2.92 | 2.52 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 1.84 | 2.36 | 3.25 | 2.85 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 1.90 | 2.03 | 2.21 | 2.03 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 14.32 | 18.66 | 13.41 | 12.59 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 4.98 | 6.65 | 7.09 | 6.46 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 18.72 | 4.24 | 6.02 | 5.36 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 3.74 | 4.31 | 4.58 | 4.44 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 21.24 | 27.99 | 21.85 | 20.30 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 7.75 | 9.91 | 11.40 | 8.97 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 5.86 | 8.32 | 10.04 | 9.17 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 5.68 | 6.74 | 7.45 | 6.94 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 28.15 | 37.33 | 30.68 | 28.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 10.74 | 14.32 | 13.55 | 13.24 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 8.58 | 11.82 | 13.39 | 12.16 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 7.56 | 9.10 | 10.46 | 9.57 |
| `bench_horner_cxx20::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 9.10 | 8.09 | 7.20 | 6.91 |
| `bench_horner_cxx20::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 36.21 | 27.32 | 17.34 | 17.94 |
| `bench_horner_cxx20::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 119.64 | 132.23 | 61.66 | 70.62 |
| `bench_horner_cxx20::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 266.60 | 332.57 | 167.88 | 180.97 |
| `bench_horner_cxx20::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 38.10 | 35.66 | 34.63 | 36.56 |
| `bench_horner_cxx20::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 296.92 | 232.56 | 205.68 | 223.54 |
| `bench_horner_cxx20::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2442.28 | 2134.97 | 1797.95 | 1885.47 |
| `bench_horner_cxx20::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 151.79 | 136.86 | 127.79 | 126.80 |
| `bench_horner_cxx20::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 2352.92 | 2022.90 | 1553.02 | 1617.93 |
| `benchMany_cxx20::13 funcs (non-many)` | — | 73.91 | 59.28 | 67.89 |

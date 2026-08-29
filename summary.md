# Benchmark Summary

All values are **ns/op** (lower is better).


| Benchmark | gcc-14 | gcc-15 | llvm-20 | llvm-21 |
| --- | ---: | ---: | ---: | ---: |
| `bench1D_cxx17::float constexpr fit` | 449.54 | 332.12 | 771.92 | 415.01 |
| `bench1D_cxx17::float coefficient-count fit` | 1036.40 | 802.00 | 929.51 | 639.31 |
| `bench1D_cxx17::float eps fit` | 2570.78 | 1720.56 | 2446.53 | 1405.31 |
| `bench1D_cxx17::float eval` | 1.95 | 1.41 | 1.99 | 1.25 |
| `bench1D_cxx17::float eval_many` | 0.26 | 0.25 | 0.17 | 0.11 |
| `bench1D_cxx17::double constexpr fit` | 494.23 | 364.93 | 805.95 | 441.36 |
| `bench1D_cxx17::double coefficient-count fit` | 1186.26 | 914.37 | 1042.21 | 681.28 |
| `bench1D_cxx17::double eps fit` | 2733.46 | 1757.73 | 2456.89 | 1385.95 |
| `bench1D_cxx17::double eval` | 2.00 | 1.41 | 1.99 | 1.25 |
| `bench1D_cxx17::double eval_many` | 0.50 | 0.50 | 0.34 | 0.21 |
| `bench1D_cxx17::complex<double> constexpr fit` | 580.91 | 464.20 | 684.99 | 402.28 |
| `bench1D_cxx17::complex<double> coefficient-count fit` | 1352.64 | 1047.60 | 1415.72 | 895.04 |
| `bench1D_cxx17::complex<double> eps fit` | 2955.98 | 1998.61 | 2749.98 | 1641.20 |
| `bench1D_cxx17::complex<double> eval` | 2.53 | 1.90 | 2.21 | 1.29 |
| `bench1D_cxx17::complex<double> eval_many` | 1.35 | 0.83 | 1.26 | 0.85 |
| `bench1D_cxx17::complex<float> constexpr fit` | 622.77 | 450.06 | 925.10 | 593.06 |
| `bench1D_cxx17::complex<float> coefficient-count fit` | 1413.47 | 1084.73 | 1649.85 | 1070.35 |
| `bench1D_cxx17::complex<float> eps fit` | 9282.11 | 7769.98 | 2836.42 | 1710.13 |
| `bench1D_cxx17::complex<float> eval` | 3.53 | 2.80 | 3.95 | 1.94 |
| `bench1D_cxx17::complex<float> eval_many` | 0.70 | 0.46 | 0.67 | 0.45 |
| `bench1D_cxx20::float constexpr fit` | 443.33 | 332.36 | 774.09 | 415.43 |
| `bench1D_cxx20::float coefficient-count fit` | 1037.15 | 793.37 | 927.83 | 634.55 |
| `bench1D_cxx20::float eps fit` | 2619.04 | 1711.72 | 2451.26 | 1412.44 |
| `bench1D_cxx20::float eval` | 2.04 | 1.41 | 1.99 | 1.25 |
| `bench1D_cxx20::float eval_many` | 0.25 | 0.27 | 0.17 | 0.11 |
| `bench1D_cxx20::double constexpr fit` | 480.03 | 354.30 | 789.29 | 444.06 |
| `bench1D_cxx20::double coefficient-count fit` | 1192.52 | 910.99 | 1047.41 | 677.12 |
| `bench1D_cxx20::double eps fit` | 2803.87 | 1765.02 | 2493.38 | 1421.55 |
| `bench1D_cxx20::double eval` | 1.96 | 1.41 | 2.01 | 1.25 |
| `bench1D_cxx20::double eval_many` | 0.50 | 0.49 | 0.34 | 0.21 |
| `bench1D_cxx20::complex<double> constexpr fit` | 567.95 | 476.96 | 670.91 | 401.50 |
| `bench1D_cxx20::complex<double> coefficient-count fit` | 1344.79 | 1101.58 | 1434.13 | 890.24 |
| `bench1D_cxx20::complex<double> eps fit` | 2920.24 | 2007.84 | 2746.83 | 1606.77 |
| `bench1D_cxx20::complex<double> eval` | 2.49 | 1.90 | 2.16 | 1.29 |
| `bench1D_cxx20::complex<double> eval_many` | 1.34 | 0.83 | 1.26 | 0.85 |
| `bench1D_cxx20::complex<float> constexpr fit` | 642.56 | 435.81 | 926.78 | 593.34 |
| `bench1D_cxx20::complex<float> coefficient-count fit` | 1411.30 | 1019.16 | 1666.36 | 1066.98 |
| `bench1D_cxx20::complex<float> eps fit` | 9268.46 | 7719.89 | 2839.56 | 1704.66 |
| `bench1D_cxx20::complex<float> eval` | 3.42 | 2.80 | 3.97 | 1.95 |
| `bench1D_cxx20::complex<float> eval_many` | 0.69 | 0.46 | 0.67 | 0.46 |
| `benchMany_cxx17::1 funcs` | 7.52 | 6.05 | 5.25 | 3.42 |
| `benchMany_cxx17::1 funcs (non-many)` | 5.21 | 3.26 | 5.09 | 2.75 |
| `benchMany_cxx17::2 funcs` | 8.45 | 7.84 | 5.26 | 3.44 |
| `benchMany_cxx17::2 funcs (non-many)` | 9.85 | 6.58 | 10.07 | 5.52 |
| `benchMany_cxx17::3 funcs` | 8.54 | 8.03 | 5.06 | 3.41 |
| `benchMany_cxx17::3 funcs (non-many)` | 16.56 | 10.59 | 13.46 | 7.82 |
| `benchMany_cxx17::4 funcs` | 8.50 | 6.56 | 5.26 | 3.42 |
| `benchMany_cxx17::4 funcs (non-many)` | 17.80 | 13.01 | 19.30 | 10.59 |
| `benchMany_cxx17::5 funcs` | 11.07 | 7.90 | 7.88 | 5.45 |
| `benchMany_cxx17::5 funcs (non-many)` | 22.69 | 15.40 | 23.56 | 13.30 |
| `benchMany_cxx17::6 funcs` | 11.12 | 8.11 | 8.02 | 5.45 |
| `benchMany_cxx17::6 funcs (non-many)` | 30.57 | 21.82 | 28.44 | 16.58 |
| `benchMany_cxx17::7 funcs` | 12.42 | 8.95 | 8.15 | 7.34 |
| `benchMany_cxx17::7 funcs (non-many)` | 32.67 | 21.69 | 32.96 | 18.38 |
| `benchMany_cxx17::8 funcs` | 11.72 | 8.05 | 8.04 | 7.11 |
| `benchMany_cxx17::8 funcs (non-many)` | 35.04 | 24.46 | 36.36 | 20.61 |
| `benchMany_cxx17::9 funcs` | 15.29 | 9.16 | 12.76 | 7.79 |
| `benchMany_cxx17::9 funcs (non-many)` | 46.19 | 31.92 | 40.58 | 22.91 |
| `benchMany_cxx17::10 funcs` | 15.32 | 9.25 | 12.71 | 20.02 |
| `benchMany_cxx17::10 funcs (non-many)` | 45.80 | 30.92 | 44.59 | 25.34 |
| `benchMany_cxx17::11 funcs` | 16.37 | 10.08 | 11.25 | 7.89 |
| `benchMany_cxx17::11 funcs (non-many)` | 57.85 | 38.53 | 48.63 | 27.80 |
| `benchMany_cxx17::12 funcs` | 15.65 | 9.31 | 13.24 | 7.66 |
| `benchMany_cxx17::12 funcs (non-many)` | 57.30 | 42.17 | 52.87 | 30.11 |
| `benchMany_cxx17::13 funcs` | 17.61 | 11.57 | — | 9.88 |
| `benchMany_cxx17::13 funcs (non-many)` | 60.55 | 40.87 | 57.09 | 33.11 |
| `benchMany_cxx17::14 funcs` | 17.67 | 11.63 | 21.21 | 9.50 |
| `benchMany_cxx17::14 funcs (non-many)` | 66.77 | 49.61 | 61.32 | 36.28 |
| `benchMany_cxx17::15 funcs` | 19.10 | 11.61 | 21.16 | 11.52 |
| `benchMany_cxx17::15 funcs (non-many)` | 74.19 | 52.20 | 68.98 | 42.69 |
| `benchMany_cxx17::16 funcs` | 17.25 | 11.53 | 14.45 | 11.81 |
| `benchMany_cxx17::16 funcs (non-many)` | 76.38 | 51.12 | 78.62 | 46.23 |
| `benchMany_cxx17::8 funcs (mixed domains)` | 11.97 | 8.01 | 8.09 | 6.26 |
| `benchMany_cxx17::8 funcs (mixed domains, non-many)` | 35.16 | 24.53 | 35.72 | 23.01 |
| `benchMany_cxx17::16 funcs (mixed domains)` | 17.21 | 11.45 | 14.76 | 11.45 |
| `benchMany_cxx17::16 funcs (mixed domains, non-many)` | 76.71 | 50.91 | 69.84 | 45.15 |
| `benchMany_cxx17::8 funcs (complex outputs, bulk)` | 38.48 | 22.44 | 25.77 | 22.25 |
| `benchMany_cxx20::1 funcs` | 7.53 | 6.14 | 5.25 | 3.90 |
| `benchMany_cxx20::1 funcs (non-many)` | 5.19 | 3.29 | 5.11 | 3.08 |
| `benchMany_cxx20::2 funcs` | 8.48 | 6.10 | 5.27 | 3.88 |
| `benchMany_cxx20::2 funcs (non-many)` | 9.83 | 6.53 | 10.06 | 6.01 |
| `benchMany_cxx20::3 funcs` | 8.55 | 8.47 | 5.07 | 3.90 |
| `benchMany_cxx20::3 funcs (non-many)` | 16.67 | 10.58 | 13.87 | 8.75 |
| `benchMany_cxx20::4 funcs` | 9.04 | 6.16 | 5.25 | 3.91 |
| `benchMany_cxx20::4 funcs (non-many)` | 17.85 | 13.00 | 19.41 | 10.67 |
| `benchMany_cxx20::5 funcs` | 10.90 | 8.42 | 8.20 | 6.24 |
| `benchMany_cxx20::5 funcs (non-many)` | 22.85 | 15.51 | 23.78 | 15.06 |
| `benchMany_cxx20::6 funcs` | 10.89 | 8.03 | 8.43 | 6.23 |
| `benchMany_cxx20::6 funcs (non-many)` | 31.21 | 21.90 | 28.43 | 18.54 |
| `benchMany_cxx20::7 funcs` | 12.06 | 8.96 | 8.06 | 6.02 |
| `benchMany_cxx20::7 funcs (non-many)` | 32.77 | 21.68 | 32.74 | 19.88 |
| `benchMany_cxx20::8 funcs` | 11.40 | 8.17 | 8.02 | 5.77 |
| `benchMany_cxx20::8 funcs (non-many)` | 36.29 | 25.19 | 35.70 | 22.98 |
| `benchMany_cxx20::9 funcs` | 14.82 | 9.26 | 12.70 | 8.22 |
| `benchMany_cxx20::9 funcs (non-many)` | 48.27 | 31.99 | 40.54 | 25.43 |
| `benchMany_cxx20::10 funcs` | 14.97 | 9.37 | 13.48 | 8.85 |
| `benchMany_cxx20::10 funcs (non-many)` | 46.59 | 31.04 | 44.28 | 25.46 |
| `benchMany_cxx20::11 funcs` | 16.24 | 10.38 | 11.98 | 7.72 |
| `benchMany_cxx20::11 funcs (non-many)` | 60.24 | 38.87 | 50.43 | 27.65 |
| `benchMany_cxx20::12 funcs` | 15.20 | 9.33 | 12.81 | 7.62 |
| `benchMany_cxx20::12 funcs (non-many)` | 59.39 | 42.56 | 54.70 | 30.08 |
| `benchMany_cxx20::13 funcs` | 17.08 | 11.56 | 14.24 | 15.86 |
| `benchMany_cxx20::13 funcs (non-many)` | 61.32 | 40.95 | 58.95 | 32.56 |
| `benchMany_cxx20::14 funcs` | 18.00 | 11.54 | 14.64 | 9.98 |
| `benchMany_cxx20::14 funcs (non-many)` | 68.10 | 49.43 | 62.88 | 34.96 |
| `benchMany_cxx20::15 funcs` | 19.19 | 11.73 | 14.71 | 18.26 |
| `benchMany_cxx20::15 funcs (non-many)` | 77.83 | 52.32 | 65.66 | 37.57 |
| `benchMany_cxx20::16 funcs` | 17.03 | 11.48 | 14.72 | 9.75 |
| `benchMany_cxx20::16 funcs (non-many)` | 77.40 | 51.14 | 69.89 | 40.29 |
| `benchMany_cxx20::8 funcs (mixed domains)` | 11.41 | 8.21 | 8.37 | 5.44 |
| `benchMany_cxx20::8 funcs (mixed domains, non-many)` | 35.35 | 25.19 | 35.72 | 20.64 |
| `benchMany_cxx20::16 funcs (mixed domains)` | 17.03 | 11.38 | 14.74 | 10.00 |
| `benchMany_cxx20::16 funcs (mixed domains, non-many)` | 76.41 | 51.67 | 69.77 | 40.35 |
| `benchMany_cxx20::8 funcs (complex outputs, bulk)` | 38.66 | 25.28 | 25.78 | 20.15 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 build` | 60713.23 | 43272.53 | 82266.87 | 40591.84 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16 build` | 88811.08 | 63358.38 | 107086.68 | 57493.32 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8 build` | 239969.09 | 167051.48 | 414870.35 | 187727.03 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8 build` | 317952.06 | 222211.19 | 528916.55 | 248788.22 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16 build` | 1520596.00 | 1099055.10 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8 build` | 2538369.25 | 1814183.33 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8 build` | 3338894.00 | 2389866.75 | — | — |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16` | 86.80 | 72.45 | 82.64 | 41.95 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16` | 131.23 | 107.87 | 171.85 | 89.65 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8` | 272.79 | 215.22 | 364.39 | 205.40 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8` | 513.98 | 286.83 | 197.15 | 138.89 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16` | 2191.49 | 2187.01 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8` | 3498.97 | 3482.88 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8` | 4764.17 | 4738.66 | — | — |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic point` | 86.84 | 72.49 | 80.25 | 42.23 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical point` | 86.88 | 72.44 | 80.07 | 42.30 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 array` | 88.35 | 72.52 | 79.82 | 41.95 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 variadic` | 86.97 | 72.42 | 81.79 | 41.96 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 manual loop` | 89660.15 | 74001.20 | 85313.52 | 44185.06 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical batch` | 36418.80 | 25065.44 | 35911.84 | 20468.28 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic loop` | 88830.04 | 74108.24 | 85235.34 | 44246.38 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 container batch` | 89652.65 | 73886.98 | 80643.43 | 42688.64 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 build` | 60874.69 | 43929.35 | 81874.60 | 40679.92 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16 build` | 88558.39 | 64575.90 | 107042.67 | 57289.35 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8 build` | 233122.07 | 166608.15 | 391787.65 | 186990.33 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8 build` | 307729.97 | 222615.84 | 527535.33 | 248314.51 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16 build` | 1523454.62 | 1120178.20 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8 build` | 2491634.50 | 1819895.50 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8 build` | 3273647.00 | 2404426.75 | — | — |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16` | 86.77 | 72.38 | 83.53 | 42.27 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16` | 130.88 | 107.80 | 178.56 | 89.75 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8` | 356.46 | 215.08 | 364.63 | 200.91 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8` | 513.75 | 286.60 | 208.33 | 127.58 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16` | 2193.04 | 2183.77 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8` | 3502.48 | 3483.09 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8` | 4780.52 | 4726.26 | — | — |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic point` | 86.80 | 72.42 | 85.39 | 41.97 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical point` | 86.84 | 72.37 | 83.75 | 41.99 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 array` | 88.36 | 72.40 | 85.05 | 42.23 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 variadic` | 86.76 | 72.37 | 83.26 | 41.53 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 manual loop` | 89965.82 | 73853.99 | 87821.41 | 43725.61 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical batch` | 36152.36 | 24483.77 | 35985.30 | 20541.91 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 span batch` | 36433.41 | 24888.31 | 35991.76 | 20538.28 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic loop` | 88835.84 | 73845.77 | 87704.24 | 43689.39 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 container batch` | 89793.94 | 73878.23 | 83731.54 | 41959.52 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.19 | 1.70 | 3.52 | 2.44 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, hybrid-runtime` | 3.13 | 2.17 | 3.17 | 1.85 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, hybrid-ct` | 2.03 | 1.52 | 1.93 | 1.20 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, hybrid-transposed` | 1.45 | 1.14 | 1.20 | 0.76 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.35 | 0.25 | 0.34 | 0.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.35 | 0.25 | 0.34 | 0.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, scalar-runtime` | 3.86 | 3.63 | 3.45 | 2.64 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, hybrid-runtime` | 3.86 | 2.51 | 3.82 | 2.24 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, hybrid-ct` | 2.97 | 2.10 | 2.47 | 1.74 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, hybrid-transposed` | 2.20 | 1.60 | 1.93 | 1.27 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, simd-aligned` | 0.49 | 0.39 | 0.52 | 0.36 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, simd-unaligned` | 0.49 | 0.40 | 0.52 | 0.37 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 5.75 | 3.75 | 4.64 | 3.80 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, hybrid-runtime` | 4.77 | 3.08 | 4.73 | 2.86 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, hybrid-ct` | 3.74 | 2.66 | 3.16 | 2.35 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, hybrid-transposed` | 2.79 | 2.16 | 2.42 | 1.54 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.66 | 0.56 | 0.71 | 0.48 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.66 | 0.56 | 0.71 | 0.48 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, scalar-runtime` | 7.88 | 4.57 | 5.27 | 3.75 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, hybrid-runtime` | 5.26 | 3.64 | 5.27 | 3.14 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, hybrid-ct` | 3.92 | 3.36 | 4.22 | 2.76 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, hybrid-transposed` | 2.98 | 2.19 | 2.79 | 1.80 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, simd-aligned` | 0.86 | 0.73 | 0.92 | 0.63 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, simd-unaligned` | 0.87 | 0.73 | 0.93 | 0.63 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 10.57 | 5.70 | 6.80 | 4.83 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, hybrid-runtime` | 5.85 | 3.95 | 6.44 | 3.41 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, hybrid-ct` | 4.62 | 3.91 | 5.63 | 3.24 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, hybrid-transposed` | 3.25 | 2.46 | 3.00 | 1.96 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 1.05 | 0.91 | 1.13 | 0.78 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 1.05 | 0.91 | 1.14 | 0.78 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 16.51 | 8.07 | 10.21 | 7.20 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, hybrid-runtime` | 7.08 | 5.50 | 7.30 | 4.81 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, hybrid-ct` | 7.48 | 5.18 | 7.55 | 4.33 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, hybrid-transposed` | 3.80 | 2.82 | 3.42 | 2.28 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.43 | 1.27 | 1.56 | 1.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.43 | 1.28 | 1.56 | — |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.35 | 2.13 | 3.61 | 3.24 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.66 | 2.73 | 3.87 | 3.32 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.35 | 2.48 | 3.59 | 3.21 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.25 | 2.43 | 3.48 | 2.87 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=12, SIMD=No, many-runtime` | 3.30 | 3.30 | 3.79 | 3.01 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=12, SIMD=No, many-runtime` | 3.85 | 3.88 | 4.25 | 3.13 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=12, SIMD=No, many-runtime` | 3.42 | 3.74 | 3.53 | 2.93 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=12, SIMD=No, many-runtime` | 3.31 | 3.62 | 3.39 | 3.21 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 4.45 | 4.33 | 4.89 | 3.97 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 5.31 | 4.83 | 5.32 | 4.08 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.97 | 4.48 | 4.78 | — |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.80 | 4.63 | 4.44 | 3.80 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=20, SIMD=No, many-runtime` | 5.92 | 4.65 | 5.38 | — |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=20, SIMD=No, many-runtime` | 6.98 | 5.41 | 6.17 | 4.16 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=20, SIMD=No, many-runtime` | 6.47 | 5.26 | 5.50 | 3.98 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=20, SIMD=No, many-runtime` | 6.45 | 5.15 | 5.40 | 3.83 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 7.57 | 5.33 | 6.88 | 5.18 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 9.27 | 6.03 | 7.65 | 5.41 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 8.38 | 6.03 | 7.10 | 5.12 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 8.24 | 5.66 | 6.98 | 5.03 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 12.45 | 7.37 | 10.31 | 7.35 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 14.41 | 7.84 | 11.47 | 7.58 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 12.86 | — | 10.39 | 7.50 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 12.89 | 7.82 | 10.66 | 7.36 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 9.33 | 6.16 | 6.19 | 3.88 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 3.89 | 1.83 | 4.03 | 1.56 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 2.35 | 1.51 | 3.11 | 2.00 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 2.19 | 1.60 | 2.07 | 1.34 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=12` | 14.01 | 8.95 | 9.93 | 6.09 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=12` | 5.62 | 2.81 | 6.40 | 2.56 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=12` | 3.80 | 2.30 | 4.70 | 3.16 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=12` | 2.99 | 2.14 | 3.01 | 1.86 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 17.53 | 11.68 | 13.74 | 8.48 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 8.15 | 4.30 | 8.04 | 3.39 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 4.89 | 3.04 | 6.21 | 4.08 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 4.40 | 3.28 | 4.34 | 2.85 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=20` | 23.34 | 14.55 | 17.31 | 10.52 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=20` | 10.60 | 5.44 | 11.43 | 4.58 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=20` | 6.92 | 5.17 | 8.06 | 5.13 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=20` | 5.23 | 3.41 | 5.48 | 2.89 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 26.23 | 17.28 | 21.31 | 12.14 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 12.99 | 6.50 | 14.27 | 5.81 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 7.93 | 4.48 | 9.17 | 5.82 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 6.79 | 5.03 | 7.46 | 4.22 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 37.40 | 22.97 | 30.65 | 16.75 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 16.40 | 7.80 | 15.64 | 7.26 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 11.47 | 6.94 | 12.58 | 7.57 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 8.72 | 5.72 | 9.73 | 5.05 |
| `bench_horner_cxx17::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 8.65 | 8.34 | 7.74 | 6.23 |
| `bench_horner_cxx17::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 27.83 | 25.45 | 23.86 | 17.12 |
| `bench_horner_cxx17::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 132.69 | 120.69 | 136.04 | 80.59 |
| `bench_horner_cxx17::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 337.75 | 209.53 | 387.18 | 214.14 |
| `bench_horner_cxx17::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 36.41 | 35.49 | 40.66 | 34.58 |
| `bench_horner_cxx17::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 231.78 | 215.69 | 283.58 | 250.93 |
| `bench_horner_cxx17::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2145.19 | 2029.36 | 2964.45 | 1995.46 |
| `bench_horner_cxx17::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 151.31 | 147.45 | 128.18 | 102.21 |
| `bench_horner_cxx17::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 1883.14 | 1704.12 | 1472.89 | 1090.80 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.22 | 1.69 | 3.82 | 2.41 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, hybrid-runtime` | 3.12 | 2.18 | 3.17 | 1.87 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, hybrid-ct` | 2.05 | 1.52 | 1.93 | 1.21 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, hybrid-transposed` | 1.45 | 1.13 | 1.19 | 0.76 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.34 | 0.28 | 0.34 | 0.21 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.36 | 0.29 | 0.34 | 0.21 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, scalar-runtime` | 3.87 | 3.63 | 3.64 | 2.64 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, hybrid-runtime` | 3.85 | 2.51 | 3.83 | 2.26 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, hybrid-ct` | 2.96 | 2.10 | 2.48 | 1.75 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, hybrid-transposed` | 2.20 | 1.60 | 1.93 | 1.20 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, simd-aligned` | 0.49 | 0.39 | 0.52 | 0.34 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, simd-unaligned` | 0.49 | 0.39 | 0.52 | 0.34 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 5.75 | 3.65 | 4.66 | 3.60 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, hybrid-runtime` | 4.78 | 3.08 | 4.75 | 2.67 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, hybrid-ct` | 3.76 | 2.66 | 3.16 | 2.23 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, hybrid-transposed` | 2.79 | 2.16 | 2.42 | 1.51 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.67 | 0.56 | 0.71 | 0.48 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.66 | 0.56 | 0.71 | 0.48 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, scalar-runtime` | 7.87 | 4.56 | 5.33 | 3.76 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, hybrid-runtime` | 5.25 | 3.64 | 5.26 | 3.15 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, hybrid-ct` | 3.87 | 3.36 | 4.22 | 2.76 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, hybrid-transposed` | 2.98 | 2.20 | 2.79 | 1.76 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, simd-aligned` | 0.87 | 0.73 | 0.92 | 0.63 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, simd-unaligned` | 0.87 | 0.73 | 0.93 | 0.63 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 10.57 | 5.70 | 6.73 | 4.83 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, hybrid-runtime` | 5.83 | 3.96 | 6.46 | 3.41 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, hybrid-ct` | 4.57 | 3.91 | 5.64 | 3.24 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, hybrid-transposed` | 3.25 | 2.46 | 3.00 | 1.92 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 1.06 | 0.91 | 1.14 | 0.78 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 1.05 | 0.91 | 1.14 | 0.78 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 16.49 | 8.06 | 10.20 | 7.19 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, hybrid-runtime` | 7.07 | 5.54 | 7.31 | 4.81 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, hybrid-ct` | 7.55 | 5.18 | 7.55 | 4.33 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, hybrid-transposed` | 3.81 | 2.82 | 3.42 | 2.25 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.44 | 1.28 | 1.56 | 1.08 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.43 | 1.29 | 1.56 | 1.09 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.33 | 2.14 | 3.60 | 2.87 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.57 | 2.77 | 3.96 | 2.97 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.34 | 2.51 | 3.64 | 2.83 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.25 | 2.58 | 3.57 | 2.86 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=12, SIMD=No, many-runtime` | 3.25 | 3.50 | 3.79 | 3.00 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=12, SIMD=No, many-runtime` | 3.84 | 3.85 | 3.93 | 3.13 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=12, SIMD=No, many-runtime` | 3.41 | 3.79 | 3.53 | 2.93 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=12, SIMD=No, many-runtime` | 3.31 | 3.69 | 3.39 | — |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 4.44 | 4.37 | 4.89 | — |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 5.31 | 5.02 | 5.00 | 4.09 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.97 | 4.66 | 4.74 | 3.89 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.81 | 4.85 | 4.44 | 3.80 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=20, SIMD=No, many-runtime` | 5.92 | 4.93 | 5.37 | 4.02 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=20, SIMD=No, many-runtime` | 6.88 | 5.61 | 5.91 | 4.16 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=20, SIMD=No, many-runtime` | 6.51 | 5.27 | 5.51 | 4.53 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=20, SIMD=No, many-runtime` | 6.45 | 5.37 | 5.40 | 4.39 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 7.56 | 5.76 | 6.88 | — |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 9.26 | 6.59 | 7.63 | 5.44 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 8.33 | 6.18 | 7.11 | 5.27 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 8.24 | 6.44 | 6.98 | 5.08 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 12.45 | 7.98 | 10.30 | 7.36 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 13.47 | 8.55 | 11.44 | 7.58 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 12.86 | 8.41 | 10.39 | 7.67 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 12.90 | 8.70 | 10.66 | 7.35 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 9.33 | 6.15 | 6.43 | 4.26 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 4.18 | 1.77 | 2.95 | 1.59 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 2.38 | 1.54 | 3.12 | 2.00 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 2.19 | 1.55 | 2.06 | 1.35 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=12` | 14.04 | 8.95 | 9.99 | 7.11 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=12` | 5.62 | 2.69 | 4.44 | 2.55 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=12` | 3.80 | 2.29 | 4.69 | 3.16 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=12` | 2.99 | 2.10 | 2.99 | 1.87 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 17.54 | 11.68 | 13.40 | 8.72 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 8.15 | 4.02 | 7.62 | 3.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 4.89 | 3.00 | 6.21 | 4.22 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 4.41 | 3.16 | 4.35 | 2.90 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=20` | 23.33 | 14.54 | 17.47 | 10.67 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=20` | 10.52 | 5.02 | 9.74 | 4.67 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=20` | 6.79 | 3.99 | 8.06 | 5.17 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=20` | 5.04 | 3.46 | 5.48 | 2.87 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 26.22 | 17.26 | 21.29 | 12.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 12.99 | 6.03 | 12.44 | 5.78 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 7.93 | 4.50 | 9.13 | 5.78 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 6.96 | 4.92 | 7.46 | 4.13 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 37.38 | 22.95 | 30.65 | 16.57 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 17.71 | 7.74 | 15.63 | 6.95 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 21.20 | 6.03 | 12.56 | 8.10 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 8.72 | 5.63 | 9.73 | 5.00 |
| `bench_horner_cxx20::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 8.62 | 8.44 | 7.75 | 6.26 |
| `bench_horner_cxx20::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 27.84 | 25.40 | 23.90 | 17.13 |
| `bench_horner_cxx20::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 132.61 | 120.68 | 136.05 | 81.12 |
| `bench_horner_cxx20::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 337.73 | 209.46 | 387.28 | 237.25 |
| `bench_horner_cxx20::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 35.97 | 35.39 | 39.71 | 34.47 |
| `bench_horner_cxx20::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 229.98 | 215.66 | 283.27 | 249.34 |
| `bench_horner_cxx20::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2097.51 | 2028.69 | 2964.19 | 1994.63 |
| `bench_horner_cxx20::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 151.21 | 160.63 | 128.21 | 102.21 |
| `bench_horner_cxx20::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 1881.09 | 1682.59 | 1469.63 | 1097.46 |
| `bench_horner_cxx17:::wavy_dash: `Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` (Unstable with ~234,512.7 iters. Increase `minEpochIterations` to e.g. 2345127)` | — | 8.31 | — | — |
| `benchMany_cxx17:::wavy_dash: `13 funcs` (Unstable with ~1,082.5 iters. Increase `minEpochIterations` to e.g. 10825)` | — | — | 20.66 | — |
| `bench_horner_cxx17:::wavy_dash: `Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` (Unstable with ~4,578.7 iters. Increase `minEpochIterations` to e.g. 45787)` | — | — | — | 1.19 |
| `bench_horner_cxx17:::wavy_dash: `Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` (Unstable with ~451,059.6 iters. Increase `minEpochIterations` to e.g. 4510596)` | — | — | — | 4.33 |
| `bench_horner_cxx17:::wavy_dash: `Dim=1, horner_many M=4, nCoeffs=20, SIMD=No, many-runtime` (Unstable with ~1,318,370.1 iters. Increase `minEpochIterations` to e.g. 13183701)` | — | — | — | 4.24 |
| `bench_horner_cxx20:::wavy_dash: `Dim=1, horner_many M=16, nCoeffs=12, SIMD=No, many-runtime` (Unstable with ~460,568.7 iters. Increase `minEpochIterations` to e.g. 4605687)` | — | — | — | 3.13 |
| `bench_horner_cxx20:::wavy_dash: `Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` (Unstable with ~1,177,702.1 iters. Increase `minEpochIterations` to e.g. 11777021)` | — | — | — | 4.49 |
| `bench_horner_cxx20:::wavy_dash: `Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` (Unstable with ~943,957.7 iters. Increase `minEpochIterations` to e.g. 9439577)` | — | — | — | 5.72 |

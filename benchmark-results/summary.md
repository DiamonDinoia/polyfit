# Benchmark Summary

All values are **ns/op** (lower is better).


| Benchmark | gcc-14 | gcc-15 | llvm-20 | llvm-21 |
| --- | ---: | ---: | ---: | ---: |
| `bench1D_cxx17::float constexpr fit` | 360.99 | 449.45 | 456.10 | 598.62 |
| `bench1D_cxx17::float coefficient-count fit` | 878.89 | 1059.36 | 691.27 | 766.39 |
| `bench1D_cxx17::float eps fit` | 2005.85 | 2528.30 | 1527.29 | 1942.44 |
| `bench1D_cxx17::float eval` | 1.50 | 2.05 | 1.35 | 1.58 |
| `bench1D_cxx17::float eval_many` | 0.18 | 0.23 | 0.11 | 0.13 |
| `bench1D_cxx17::double constexpr fit` | 405.79 | 515.82 | 473.31 | 610.31 |
| `bench1D_cxx17::double coefficient-count fit` | 957.46 | 1202.24 | 725.98 | 825.55 |
| `bench1D_cxx17::double eps fit` | 2088.59 | 2659.96 | 1523.00 | 1946.91 |
| `bench1D_cxx17::double eval` | 1.47 | 2.03 | 1.35 | 1.56 |
| `bench1D_cxx17::double eval_many` | 0.42 | 0.46 | 0.23 | 0.27 |
| `bench1D_cxx17::complex<double> constexpr fit` | 489.69 | 618.75 | 435.04 | 614.66 |
| `bench1D_cxx17::complex<double> coefficient-count fit` | 1108.65 | 1399.78 | 962.04 | 1122.04 |
| `bench1D_cxx17::complex<double> eps fit` | 2273.63 | 3004.73 | 1754.81 | 2187.40 |
| `bench1D_cxx17::complex<double> eval` | 2.07 | 2.55 | 1.38 | 1.77 |
| `bench1D_cxx17::complex<double> eval_many` | 0.95 | 1.23 | 0.86 | 0.98 |
| `bench1D_cxx17::complex<float> constexpr fit` | 502.84 | 595.71 | 639.19 | 725.08 |
| `bench1D_cxx17::complex<float> coefficient-count fit` | 1127.66 | 1398.32 | 1162.17 | 1306.25 |
| `bench1D_cxx17::complex<float> eps fit` | 7977.72 | 9664.53 | 1863.06 | 2239.94 |
| `bench1D_cxx17::complex<float> eval` | 2.75 | 3.31 | 2.13 | 3.07 |
| `bench1D_cxx17::complex<float> eval_many` | 0.53 | 0.69 | 0.47 | 0.52 |
| `bench1D_cxx20::float constexpr fit` | 351.62 | 448.52 | 454.18 | 604.25 |
| `bench1D_cxx20::float coefficient-count fit` | 842.44 | 1051.77 | 691.81 | 771.12 |
| `bench1D_cxx20::float eps fit` | 1985.97 | 2561.53 | 1522.06 | 1941.28 |
| `bench1D_cxx20::float eval` | 1.49 | 2.04 | 1.35 | 1.55 |
| `bench1D_cxx20::float eval_many` | 0.18 | 0.23 | 0.11 | 0.13 |
| `bench1D_cxx20::double constexpr fit` | 391.18 | 504.36 | 473.35 | 608.45 |
| `bench1D_cxx20::double coefficient-count fit` | 952.50 | 1190.25 | 741.06 | 827.71 |
| `bench1D_cxx20::double eps fit` | 2069.77 | 2676.46 | 1493.81 | 1938.69 |
| `bench1D_cxx20::double eval` | 1.47 | 2.08 | 1.35 | 1.57 |
| `bench1D_cxx20::double eval_many` | 0.36 | 0.46 | 0.23 | 0.26 |
| `bench1D_cxx20::complex<double> constexpr fit` | 475.77 | 607.87 | 438.18 | 525.58 |
| `bench1D_cxx20::complex<double> coefficient-count fit` | 1103.65 | 1457.70 | 976.29 | 1098.23 |
| `bench1D_cxx20::complex<double> eps fit` | 2257.17 | 3018.56 | 1758.86 | 2206.73 |
| `bench1D_cxx20::complex<double> eval` | 2.08 | 2.65 | 1.38 | 1.73 |
| `bench1D_cxx20::complex<double> eval_many` | 0.95 | 1.26 | 0.86 | 0.92 |
| `bench1D_cxx20::complex<float> constexpr fit` | 509.00 | 568.62 | 637.44 | 724.78 |
| `bench1D_cxx20::complex<float> coefficient-count fit` | 1122.36 | 1361.20 | 1169.55 | 1317.31 |
| `bench1D_cxx20::complex<float> eps fit` | 7958.08 | 9582.65 | 1871.57 | 2244.52 |
| `bench1D_cxx20::complex<float> eval` | 2.77 | 3.33 | 2.13 | 3.06 |
| `bench1D_cxx20::complex<float> eval_many` | 0.53 | 0.69 | 0.45 | 0.47 |
| `benchMany_cxx17::1 funcs` | 5.75 | 7.11 | 3.67 | 4.08 |
| `benchMany_cxx17::1 funcs (non-many)` | 3.90 | 5.15 | 2.80 | 3.94 |
| `benchMany_cxx17::2 funcs` | 5.92 | 7.40 | 3.70 | 4.09 |
| `benchMany_cxx17::2 funcs (non-many)` | 7.62 | 9.83 | 5.99 | 7.80 |
| `benchMany_cxx17::3 funcs` | 6.08 | 7.72 | 3.68 | 3.94 |
| `benchMany_cxx17::3 funcs (non-many)` | 12.57 | 15.92 | 8.03 | 10.44 |
| `benchMany_cxx17::4 funcs` | 6.20 | 7.96 | 3.69 | 4.11 |
| `benchMany_cxx17::4 funcs (non-many)` | 14.70 | 19.24 | 10.91 | 14.42 |
| `benchMany_cxx17::5 funcs` | 9.53 | 12.37 | 5.59 | 6.24 |
| `benchMany_cxx17::5 funcs (non-many)` | 17.11 | 22.24 | 14.08 | 18.31 |
| `benchMany_cxx17::6 funcs` | 9.46 | 12.24 | 5.61 | 6.40 |
| `benchMany_cxx17::6 funcs (non-many)` | 25.07 | 32.72 | 16.60 | 22.05 |
| `benchMany_cxx17::7 funcs` | 9.51 | 13.02 | 5.60 | 6.22 |
| `benchMany_cxx17::7 funcs (non-many)` | 25.12 | 32.31 | 18.83 | 25.47 |
| `benchMany_cxx17::8 funcs` | 9.76 | 13.25 | 5.61 | 6.14 |
| `benchMany_cxx17::8 funcs (non-many)` | 28.05 | 36.19 | 21.20 | 27.72 |
| `benchMany_cxx17::9 funcs` | 13.79 | 15.83 | 8.47 | 10.32 |
| `benchMany_cxx17::9 funcs (non-many)` | 35.07 | 45.84 | 23.57 | 31.30 |
| `benchMany_cxx17::10 funcs` | 26.17 | 16.31 | 8.50 | 10.60 |
| `benchMany_cxx17::10 funcs (non-many)` | 34.61 | 44.95 | 26.04 | 34.40 |
| `benchMany_cxx17::11 funcs` | 12.52 | 16.62 | 8.44 | 9.23 |
| `benchMany_cxx17::11 funcs (non-many)` | 43.16 | 56.98 | 28.79 | 37.74 |
| `benchMany_cxx17::12 funcs` | 19.17 | 16.02 | 8.17 | 10.22 |
| `benchMany_cxx17::12 funcs (non-many)` | 47.00 | 62.80 | 30.90 | 40.98 |
| `benchMany_cxx17::13 funcs` | 15.04 | 19.30 | 10.79 | 11.03 |
| `benchMany_cxx17::13 funcs (non-many)` | 46.47 | 55.83 | 33.93 | 44.43 |
| `benchMany_cxx17::14 funcs` | 19.94 | 18.99 | 10.55 | 11.39 |
| `benchMany_cxx17::14 funcs (non-many)` | 54.44 | 72.31 | 36.52 | 47.55 |
| `benchMany_cxx17::15 funcs` | 14.88 | 19.12 | 11.31 | 11.41 |
| `benchMany_cxx17::15 funcs (non-many)` | 57.87 | 76.02 | 38.54 | 50.87 |
| `benchMany_cxx17::16 funcs` | 14.52 | 19.10 | 10.86 | 11.48 |
| `benchMany_cxx17::16 funcs (non-many)` | 53.22 | 68.62 | 41.10 | 54.20 |
| `benchMany_cxx17::8 funcs (mixed domains)` | 9.44 | 12.20 | 5.59 | 6.23 |
| `benchMany_cxx17::8 funcs (mixed domains, non-many)` | 28.07 | 36.45 | 21.21 | 28.62 |
| `benchMany_cxx17::16 funcs (mixed domains)` | 14.46 | 18.66 | 11.14 | 11.32 |
| `benchMany_cxx17::16 funcs (mixed domains, non-many)` | 53.16 | 68.72 | 41.15 | 55.43 |
| `benchMany_cxx17::8 funcs (complex outputs, bulk)` | 30.72 | 39.91 | 21.59 | 22.19 |
| `benchMany_cxx20::1 funcs` | 5.75 | 7.10 | 3.68 | 4.10 |
| `benchMany_cxx20::1 funcs (non-many)` | 3.94 | 5.10 | 2.85 | 3.96 |
| `benchMany_cxx20::2 funcs` | 5.96 | 7.76 | 3.68 | 4.09 |
| `benchMany_cxx20::2 funcs (non-many)` | 7.60 | 9.87 | 6.07 | 7.83 |
| `benchMany_cxx20::3 funcs` | 6.06 | 7.82 | 3.70 | 3.93 |
| `benchMany_cxx20::3 funcs (non-many)` | 12.56 | 15.99 | 8.64 | 10.80 |
| `benchMany_cxx20::4 funcs` | 6.05 | 7.65 | 3.68 | 4.12 |
| `benchMany_cxx20::4 funcs (non-many)` | 14.71 | 19.39 | 11.93 | 15.09 |
| `benchMany_cxx20::5 funcs` | 9.52 | 12.40 | 5.61 | 6.21 |
| `benchMany_cxx20::5 funcs (non-many)` | 17.09 | 22.26 | 14.22 | 18.46 |
| `benchMany_cxx20::6 funcs` | 9.35 | 12.17 | 5.59 | 6.21 |
| `benchMany_cxx20::6 funcs (non-many)` | 25.46 | 32.72 | 17.22 | 22.04 |
| `benchMany_cxx20::7 funcs` | 9.51 | 12.42 | 5.55 | 6.21 |
| `benchMany_cxx20::7 funcs (non-many)` | 25.11 | 32.26 | 19.17 | 25.46 |
| `benchMany_cxx20::8 funcs` | 9.40 | 12.71 | 5.59 | 6.17 |
| `benchMany_cxx20::8 funcs (non-many)` | 28.00 | 36.13 | 21.19 | 27.75 |
| `benchMany_cxx20::9 funcs` | 12.07 | 14.92 | 8.36 | 10.31 |
| `benchMany_cxx20::9 funcs (non-many)` | 35.06 | 45.90 | 23.60 | 31.18 |
| `benchMany_cxx20::10 funcs` | 12.03 | 15.36 | 8.46 | 10.54 |
| `benchMany_cxx20::10 funcs (non-many)` | 35.04 | 45.69 | 26.08 | 34.38 |
| `benchMany_cxx20::11 funcs` | 12.82 | 16.26 | 8.18 | 9.19 |
| `benchMany_cxx20::11 funcs (non-many)` | 43.28 | 57.07 | 28.48 | 37.64 |
| `benchMany_cxx20::12 funcs` | 12.36 | 17.41 | 8.13 | 10.18 |
| `benchMany_cxx20::12 funcs (non-many)` | 47.52 | 62.16 | 30.93 | 41.08 |
| `benchMany_cxx20::13 funcs` | 14.66 | 42.77 | 10.80 | 11.54 |
| `benchMany_cxx20::13 funcs (non-many)` | 46.47 | 59.92 | 33.76 | 44.31 |
| `benchMany_cxx20::14 funcs` | 15.18 | 19.27 | 10.86 | 12.16 |
| `benchMany_cxx20::14 funcs (non-many)` | 54.79 | 71.66 | 35.94 | 47.62 |
| `benchMany_cxx20::15 funcs` | 14.70 | 19.19 | 10.37 | 12.01 |
| `benchMany_cxx20::15 funcs (non-many)` | 58.14 | 75.58 | 38.72 | 51.00 |
| `benchMany_cxx20::16 funcs` | 14.47 | 18.75 | 10.78 | 11.08 |
| `benchMany_cxx20::16 funcs (non-many)` | 53.18 | 73.65 | 41.31 | 54.26 |
| `benchMany_cxx20::8 funcs (mixed domains)` | 9.34 | 23.08 | 5.59 | 6.40 |
| `benchMany_cxx20::8 funcs (mixed domains, non-many)` | 28.11 | 36.15 | 21.21 | 27.81 |
| `benchMany_cxx20::16 funcs (mixed domains)` | 14.53 | 18.52 | 10.71 | 11.11 |
| `benchMany_cxx20::16 funcs (mixed domains, non-many)` | 53.22 | 73.60 | 41.28 | 54.25 |
| `benchMany_cxx20::8 funcs (complex outputs, bulk)` | 30.78 | 39.47 | 21.63 | 21.05 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 build` | 49163.57 | 62993.24 | 43787.85 | 63835.38 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16 build` | 71439.56 | 91119.06 | 61643.26 | 82848.34 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8 build` | 194118.87 | 246568.74 | 201248.93 | 320787.28 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8 build` | 256321.13 | 327051.69 | 267490.00 | 406268.80 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16 build` | 1215304.78 | 1563912.29 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8 build` | 2033878.00 | 2690133.00 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8 build` | 2688972.75 | 3551752.33 | — | — |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16` | 77.42 | 99.20 | 43.44 | 62.04 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16` | 113.89 | 149.17 | 103.00 | 133.69 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8` | 239.55 | 307.73 | 216.60 | 282.97 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8` | 320.91 | 415.11 | 131.61 | 157.35 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16` | 1692.72 | 2199.16 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8` | 2692.19 | 3497.50 | — | — |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8` | 3657.36 | 4746.38 | — | — |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic point` | 77.61 | 98.86 | 43.41 | 62.43 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical point` | 77.55 | 98.68 | 43.39 | 64.13 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 array` | 77.60 | 99.17 | 43.14 | 61.95 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 variadic` | 75.65 | 98.82 | 43.13 | 61.72 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 manual loop` | 77712.74 | 101350.60 | 45482.89 | 66138.39 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical batch` | 28789.98 | 35526.64 | 24004.44 | 29146.41 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic loop` | 77556.97 | 101189.00 | 45536.25 | 66223.12 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 container batch` | 77979.41 | 101400.22 | 43847.24 | 62680.94 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 build` | 49115.60 | 65097.23 | 43698.82 | 64453.73 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16 build` | 71794.94 | 92097.41 | 61788.43 | 89508.93 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8 build` | 190810.83 | 267524.61 | 201979.47 | 321159.00 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8 build` | 253619.11 | 354812.66 | 268056.82 | 405564.54 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16 build` | 1234221.12 | 1627334.00 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8 build` | 2001043.00 | 2729718.00 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8 build` | 2636570.00 | 3606709.00 | — | — |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16` | 77.47 | 97.55 | 43.22 | 64.07 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16` | 114.18 | 149.07 | 96.95 | 133.54 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8` | 240.62 | 307.77 | 218.03 | 283.10 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8` | 321.52 | 415.80 | 131.47 | 169.77 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16` | 1690.09 | 2175.00 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8` | 2714.56 | 3466.16 | — | — |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8` | 3645.01 | 4705.23 | — | — |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic point` | 77.34 | 97.54 | 43.52 | 62.18 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical point` | 77.65 | 97.66 | 43.45 | 62.29 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 array` | 77.73 | 97.41 | 43.19 | 62.18 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 variadic` | 75.45 | 97.43 | 43.17 | 61.76 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 manual loop` | 77549.01 | 100028.83 | 45577.72 | 66171.99 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical batch` | 28794.22 | 35578.67 | 22810.96 | 29744.91 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 span batch` | 28794.35 | 35513.20 | 23085.22 | 29619.52 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic loop` | 77605.77 | 100034.90 | 45550.35 | 66092.39 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 container batch` | 78004.72 | 100515.46 | 43996.48 | 62470.58 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 1.91 | 2.40 | 2.60 | 2.73 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, hybrid-runtime` | 2.32 | 2.99 | 2.00 | 2.46 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, hybrid-ct` | 1.66 | 2.17 | 1.30 | 1.50 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, hybrid-transposed` | 1.06 | 1.35 | 0.82 | 0.94 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.27 | 0.34 | 0.23 | 0.27 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.27 | 0.36 | 0.23 | 0.28 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, scalar-runtime` | 3.00 | 3.21 | 2.85 | 2.90 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, hybrid-runtime` | 2.85 | 3.68 | 2.41 | 2.97 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=No, hybrid-ct` | 2.45 | 3.16 | 1.88 | 1.97 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, hybrid-transposed` | 1.64 | 2.08 | 1.31 | 1.52 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, simd-aligned` | 0.40 | 0.52 | 0.36 | 0.40 |
| `bench_horner_cxx17::Dim=1, nCoeffs=12, SIMD=Yes, simd-unaligned` | 0.40 | 0.52 | 0.37 | 0.40 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.23 | 4.23 | 3.88 | 3.58 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, hybrid-runtime` | 3.94 | 5.11 | 2.74 | 3.68 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, hybrid-ct` | 3.28 | 4.29 | 2.40 | 2.46 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, hybrid-transposed` | 2.05 | 2.59 | 1.65 | 1.92 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.55 | 0.71 | 0.52 | 0.55 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.55 | 0.71 | 0.52 | 0.55 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, scalar-runtime` | 5.57 | 5.73 | 4.25 | 4.13 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, hybrid-runtime` | 4.48 | 5.81 | 3.37 | 4.09 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=No, hybrid-ct` | 3.27 | 4.08 | 2.98 | 3.28 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, hybrid-transposed` | 2.27 | 2.94 | 1.93 | 2.15 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, simd-aligned` | 0.72 | 0.92 | 0.68 | 0.72 |
| `bench_horner_cxx17::Dim=1, nCoeffs=20, SIMD=Yes, simd-unaligned` | 0.73 | 0.93 | 0.68 | 0.72 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.22 | 7.65 | 5.49 | 5.32 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, hybrid-runtime` | 5.29 | 7.21 | 3.66 | 4.98 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, hybrid-ct` | 3.86 | 4.95 | 3.49 | 4.39 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, hybrid-transposed` | 2.48 | 3.18 | 2.07 | 2.36 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 0.88 | 1.14 | 0.85 | 0.88 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 0.90 | 1.14 | 0.85 | 0.88 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 11.17 | 11.47 | 8.89 | 7.92 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, hybrid-runtime` | 5.49 | 7.10 | 5.17 | 5.66 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, hybrid-ct` | 6.01 | 7.75 | 4.66 | 5.87 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, hybrid-transposed` | 2.77 | 3.57 | 2.43 | 2.66 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.21 | 1.56 | 1.18 | 1.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.21 | 1.56 | 1.18 | 1.21 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.11 | 2.45 | 3.01 | 2.96 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.12 | 2.51 | 3.38 | 3.11 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.06 | 2.42 | 3.17 | 2.84 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 1.81 | 2.16 | 3.18 | 2.77 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=12, SIMD=No, many-runtime` | 2.58 | 3.11 | 3.11 | 3.11 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=12, SIMD=No, many-runtime` | 2.80 | 3.35 | 3.40 | 3.01 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=12, SIMD=No, many-runtime` | 2.66 | 3.32 | 3.18 | 2.74 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=12, SIMD=No, many-runtime` | 2.51 | 3.05 | 3.06 | 2.56 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 3.58 | 4.43 | 4.15 | 3.89 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 3.70 | 4.65 | 4.44 | 3.82 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 3.59 | 4.37 | 4.22 | 3.62 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 3.38 | 4.26 | 4.11 | 3.47 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=20, SIMD=No, many-runtime` | 4.81 | 5.87 | 4.30 | 4.40 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=20, SIMD=No, many-runtime` | 4.98 | 6.06 | 4.61 | 4.51 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=20, SIMD=No, many-runtime` | 4.74 | 5.93 | 4.36 | 4.22 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=20, SIMD=No, many-runtime` | 4.61 | 5.59 | 4.19 | 4.15 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 6.11 | 7.55 | 5.55 | 5.50 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 6.39 | 7.82 | 5.93 | 5.82 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 6.13 | 7.80 | 5.60 | 5.45 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 5.88 | 7.40 | 5.50 | 5.44 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 9.06 | 11.39 | 8.02 | 8.03 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 9.38 | 11.86 | 8.43 | 8.53 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 9.18 | 11.63 | 8.32 | 8.06 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 9.00 | 11.26 | 8.01 | 8.27 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 8.30 | 10.65 | 4.36 | 4.99 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.54 | 3.04 | 1.64 | 2.79 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 1.95 | 2.16 | 2.06 | 2.41 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 1.89 | 2.20 | 1.38 | 1.64 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=12` | 12.35 | 15.82 | 7.50 | 7.76 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=12` | 4.22 | 5.16 | 2.62 | 3.40 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=12` | 2.95 | 3.76 | 3.25 | 3.64 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=12` | 2.60 | 3.09 | 1.90 | 2.33 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 15.46 | 19.83 | 9.32 | 10.40 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 5.79 | 7.19 | 3.49 | 5.71 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 3.92 | 4.98 | 4.19 | 4.82 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 3.73 | 4.48 | 2.99 | 3.37 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=20` | 20.51 | 26.54 | 11.36 | 13.61 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=20` | 7.40 | 9.37 | 5.04 | 7.51 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=20` | 5.13 | 6.61 | 5.31 | 6.26 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=20` | 4.00 | 5.22 | 2.94 | 4.25 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 23.09 | 29.68 | 13.24 | 16.54 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 9.12 | 11.22 | 6.22 | 9.65 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 5.95 | 7.29 | 5.95 | 7.10 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 5.59 | 6.96 | 4.44 | 5.73 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 32.78 | 42.18 | 17.84 | 23.80 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 12.06 | 16.05 | 7.48 | 12.14 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 8.45 | 10.96 | 7.76 | 9.74 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 7.32 | 9.52 | 5.14 | 7.55 |
| `bench_horner_cxx17::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 7.01 | 8.18 | 6.77 | 6.01 |
| `bench_horner_cxx17::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 22.28 | 26.99 | 18.46 | 19.19 |
| `bench_horner_cxx17::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 86.13 | 112.56 | 89.14 | 105.60 |
| `bench_horner_cxx17::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 229.81 | 297.92 | 236.66 | 300.48 |
| `bench_horner_cxx17::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 30.25 | 37.95 | 36.00 | 30.82 |
| `bench_horner_cxx17::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 183.74 | 233.30 | 271.84 | 221.57 |
| `bench_horner_cxx17::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 1485.60 | 1916.20 | 2286.46 | 2316.81 |
| `bench_horner_cxx17::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 104.05 | 145.21 | 105.15 | 99.29 |
| `bench_horner_cxx17::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 1375.21 | 1846.66 | 1185.28 | 1138.32 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 1.91 | 2.38 | 2.60 | 2.96 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, hybrid-runtime` | 2.32 | 2.99 | 1.98 | 2.46 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, hybrid-ct` | 1.67 | 2.16 | 1.30 | 1.51 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, hybrid-transposed` | 1.06 | 1.35 | 0.82 | 0.95 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.26 | 0.34 | 0.23 | 0.26 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.27 | 0.34 | 0.23 | 0.27 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, scalar-runtime` | 3.00 | 3.22 | 2.85 | 2.83 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, hybrid-runtime` | 2.85 | 3.68 | 2.42 | 2.96 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=No, hybrid-ct` | 2.43 | 3.13 | 1.88 | 1.96 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, hybrid-transposed` | 1.64 | 2.08 | 1.31 | 1.54 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, simd-aligned` | 0.40 | 0.52 | 0.36 | 0.40 |
| `bench_horner_cxx20::Dim=1, nCoeffs=12, SIMD=Yes, simd-unaligned` | 0.40 | 0.52 | 0.37 | 0.40 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.20 | 4.21 | 3.89 | 3.58 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, hybrid-runtime` | 3.95 | 5.08 | 2.74 | 3.69 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, hybrid-ct` | 3.30 | 4.26 | 2.41 | 2.49 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, hybrid-transposed` | 2.05 | 2.58 | 1.65 | 1.93 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.55 | 0.71 | 0.52 | 0.55 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.55 | 0.71 | 0.52 | 0.55 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, scalar-runtime` | 5.58 | 5.72 | 4.25 | 4.17 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, hybrid-runtime` | 4.49 | 5.80 | 3.40 | 4.08 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=No, hybrid-ct` | 3.27 | 4.10 | 2.98 | 3.27 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, hybrid-transposed` | 2.27 | 2.94 | 1.93 | 2.16 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, simd-aligned` | 0.72 | 0.93 | 0.68 | 0.72 |
| `bench_horner_cxx20::Dim=1, nCoeffs=20, SIMD=Yes, simd-unaligned` | 0.72 | 0.93 | 0.68 | 0.72 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.26 | 7.66 | 5.49 | 5.35 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, hybrid-runtime` | 5.29 | 6.83 | 3.67 | 5.50 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, hybrid-ct` | 3.84 | 4.95 | 3.49 | 4.49 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, hybrid-transposed` | 2.48 | 3.18 | 2.11 | 2.36 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 0.88 | 1.14 | 0.85 | 0.88 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 0.88 | 1.14 | 0.85 | 0.88 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 11.20 | 11.51 | 8.89 | 7.92 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, hybrid-runtime` | 5.55 | 7.13 | 5.18 | 5.67 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, hybrid-ct` | 6.00 | 7.74 | 4.67 | 5.88 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, hybrid-transposed` | 2.76 | 3.57 | 2.44 | 2.67 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.21 | 1.56 | 1.18 | 1.21 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.21 | 1.56 | 1.18 | 1.21 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.12 | 2.78 | 3.02 | 3.01 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.13 | 2.65 | 3.15 | 3.02 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.02 | 2.45 | 2.97 | 2.80 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 1.81 | 2.30 | 2.93 | 2.68 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=12, SIMD=No, many-runtime` | 2.57 | 3.19 | 3.11 | 3.11 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=12, SIMD=No, many-runtime` | 2.80 | 3.63 | 3.40 | 3.01 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=12, SIMD=No, many-runtime` | 2.66 | 3.28 | 3.17 | 2.73 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=12, SIMD=No, many-runtime` | 2.50 | 3.28 | 3.06 | 2.56 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 3.57 | 4.74 | 4.14 | 3.88 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 3.70 | 4.67 | 4.44 | 3.81 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 3.59 | 4.41 | 4.21 | 3.65 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 3.37 | 4.44 | 4.10 | 3.40 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=20, SIMD=No, many-runtime` | 4.81 | 6.10 | 4.45 | 4.34 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=20, SIMD=No, many-runtime` | 4.99 | 6.23 | 4.97 | 4.58 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=20, SIMD=No, many-runtime` | 4.74 | 5.93 | 4.49 | 4.23 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=20, SIMD=No, many-runtime` | 4.61 | 5.82 | 4.35 | 4.14 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 6.10 | 7.71 | 5.74 | 5.44 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 6.38 | 7.93 | 6.21 | 5.82 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 6.32 | 7.85 | 5.73 | 5.47 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 5.87 | 7.56 | 5.73 | 5.43 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 9.04 | 11.50 | 8.83 | 8.02 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 9.37 | 12.09 | 8.90 | 8.53 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 9.17 | 11.76 | 8.89 | 8.06 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 8.99 | 11.70 | 8.73 | 8.26 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 8.29 | 10.65 | 3.99 | 4.82 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.52 | 3.04 | 1.61 | 2.79 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 1.82 | 2.16 | 2.05 | 2.42 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 1.88 | 2.41 | 1.38 | 1.60 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=12` | 12.34 | 15.82 | 6.28 | 7.70 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=12` | 4.22 | 5.16 | 2.61 | 3.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=12` | 2.94 | 3.76 | 3.26 | 3.64 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=12` | 2.59 | 3.09 | 1.91 | 2.33 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 15.48 | 19.81 | 9.15 | 10.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 5.79 | 7.19 | 3.48 | 5.70 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 3.92 | 4.92 | 4.16 | 4.82 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 3.73 | 4.94 | 2.93 | 3.37 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=20` | 20.51 | 26.37 | 11.32 | 13.42 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=20` | 7.39 | 9.36 | 4.93 | 7.51 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=20` | 5.13 | 6.61 | 5.34 | 6.33 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=20` | 3.96 | 5.23 | 2.96 | 4.25 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 23.10 | 29.69 | 13.10 | 16.53 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 9.05 | 11.23 | 6.30 | 9.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 5.96 | 7.27 | 5.97 | 7.11 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 5.59 | 7.51 | 4.37 | 5.79 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 32.79 | 42.18 | 18.10 | 23.79 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 12.07 | 16.06 | 7.66 | 12.13 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 8.44 | 11.73 | 7.75 | 9.77 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 7.15 | 9.16 | 5.19 | 7.55 |
| `bench_horner_cxx20::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 7.02 | 8.26 | 6.76 | 6.12 |
| `bench_horner_cxx20::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 22.29 | 27.04 | 18.49 | 19.21 |
| `bench_horner_cxx20::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 86.16 | 112.54 | 87.10 | 105.61 |
| `bench_horner_cxx20::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 229.72 | 298.19 | 232.40 | 300.66 |
| `bench_horner_cxx20::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 30.21 | 38.42 | 37.31 | 31.56 |
| `bench_horner_cxx20::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 183.51 | 240.93 | 271.02 | 221.75 |
| `bench_horner_cxx20::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 1481.39 | 1917.24 | 2164.43 | 2323.96 |
| `bench_horner_cxx20::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 103.83 | 145.56 | 105.21 | 99.37 |
| `bench_horner_cxx20::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 1392.34 | 1907.55 | 1165.98 | 1135.83 |

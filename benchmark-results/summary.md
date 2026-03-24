# Benchmark Summary

All values are **ns/op** (lower is better).


| Benchmark | gcc-14 | gcc-15 | llvm-20 | llvm-21 |
| --- | ---: | ---: | ---: | ---: |
| `bench1D_cxx17::float constexpr fit` | 466.03 | 457.93 | 732.55 | 772.31 |
| `bench1D_cxx17::float coefficient-count fit` | 1016.66 | 952.29 | 941.22 | 928.54 |
| `bench1D_cxx17::float eps fit` | 2355.21 | 2284.88 | 2260.92 | 2397.71 |
| `bench1D_cxx17::float eval` | 1.75 | 1.73 | 2.07 | 1.83 |
| `bench1D_cxx17::float eval_many` | 0.25 | 0.26 | 0.17 | 0.17 |
| `bench1D_cxx17::double constexpr fit` | 504.36 | 501.09 | 750.41 | 795.55 |
| `bench1D_cxx17::double coefficient-count fit` | 1193.65 | 1091.09 | 1034.32 | 1050.88 |
| `bench1D_cxx17::double eps fit` | 2457.45 | 2415.79 | 2247.33 | 2395.44 |
| `bench1D_cxx17::double eval` | 2.03 | 2.04 | 2.02 | 1.79 |
| `bench1D_cxx17::double eval_many` | 0.51 | 0.52 | 0.35 | 0.33 |
| `bench1D_cxx17::complex<double> constexpr fit` | 581.61 | 586.83 | 634.20 | 671.20 |
| `bench1D_cxx17::complex<double> coefficient-count fit` | 1388.97 | 1390.56 | 1463.64 | 1535.22 |
| `bench1D_cxx17::complex<double> eps fit` | 2669.57 | 2706.63 | 2553.20 | 2712.37 |
| `bench1D_cxx17::complex<double> eval` | 2.18 | 2.40 | 2.08 | 1.86 |
| `bench1D_cxx17::complex<double> eval_many` | 1.35 | 1.34 | 1.40 | 1.32 |
| `bench1D_cxx17::complex<float> constexpr fit` | 622.12 | 584.39 | 852.12 | 922.21 |
| `bench1D_cxx17::complex<float> coefficient-count fit` | 1311.67 | 1269.75 | 1572.98 | 1643.96 |
| `bench1D_cxx17::complex<float> eps fit` | 2660.68 | 2635.48 | 2644.42 | 2754.20 |
| `bench1D_cxx17::complex<float> eval` | 3.90 | 3.88 | 3.43 | 3.51 |
| `bench1D_cxx17::complex<float> eval_many` | 0.68 | 0.68 | 0.66 | 0.67 |
| `bench1D_cxx20::float constexpr fit` | 442.35 | 436.05 | 733.98 | 774.03 |
| `bench1D_cxx20::float coefficient-count fit` | 1027.61 | 1008.02 | 944.79 | 923.84 |
| `bench1D_cxx20::float eps fit` | 2343.71 | 2349.27 | 2253.40 | 2368.67 |
| `bench1D_cxx20::float eval` | 1.69 | 2.03 | 2.11 | 1.76 |
| `bench1D_cxx20::float eval_many` | 0.25 | 0.26 | 0.17 | 0.17 |
| `bench1D_cxx20::double constexpr fit` | 481.19 | 476.08 | 750.14 | 800.88 |
| `bench1D_cxx20::double coefficient-count fit` | 1193.50 | 1174.33 | 1028.55 | 1041.87 |
| `bench1D_cxx20::double eps fit` | 2525.02 | 2515.38 | 2255.04 | 2385.16 |
| `bench1D_cxx20::double eval` | 1.74 | 1.73 | 2.09 | 1.87 |
| `bench1D_cxx20::double eval_many` | 0.51 | 0.50 | 0.38 | 0.34 |
| `bench1D_cxx20::complex<double> constexpr fit` | 560.85 | 581.61 | 639.30 | 671.10 |
| `bench1D_cxx20::complex<double> coefficient-count fit` | 1360.43 | 1392.30 | 1458.38 | 1622.03 |
| `bench1D_cxx20::complex<double> eps fit` | 2776.85 | 2776.33 | 2562.56 | 2689.84 |
| `bench1D_cxx20::complex<double> eval` | 2.47 | 2.33 | 2.11 | 1.80 |
| `bench1D_cxx20::complex<double> eval_many` | 1.36 | 1.35 | 1.39 | 1.19 |
| `bench1D_cxx20::complex<float> constexpr fit` | 641.24 | 574.51 | 855.63 | 924.11 |
| `bench1D_cxx20::complex<float> coefficient-count fit` | 1398.50 | 1401.95 | 1572.64 | 1647.58 |
| `bench1D_cxx20::complex<float> eps fit` | 2736.24 | 2706.86 | 2610.51 | 2728.26 |
| `bench1D_cxx20::complex<float> eval` | 3.87 | 3.89 | 3.51 | 3.47 |
| `bench1D_cxx20::complex<float> eval_many` | 0.68 | 0.68 | 0.67 | 0.60 |
| `benchMany_cxx17::1 funcs` | 7.55 | 7.64 | 5.66 | 5.28 |
| `benchMany_cxx17::1 funcs (non-many)` | 4.59 | 4.37 | 5.15 | 4.58 |
| `benchMany_cxx17::2 funcs` | 8.49 | 7.97 | 5.88 | 5.28 |
| `benchMany_cxx17::2 funcs (non-many)` | 9.38 | 9.28 | 10.00 | 8.56 |
| `benchMany_cxx17::3 funcs` | 8.61 | 8.63 | 6.14 | 5.08 |
| `benchMany_cxx17::3 funcs (non-many)` | 19.28 | 13.46 | 14.82 | 12.87 |
| `benchMany_cxx17::4 funcs` | 8.47 | 8.57 | 5.66 | 5.24 |
| `benchMany_cxx17::4 funcs (non-many)` | 18.32 | 18.04 | 19.34 | 17.02 |
| `benchMany_cxx17::5 funcs` | 10.70 | 10.62 | 8.15 | 7.93 |
| `benchMany_cxx17::5 funcs (non-many)` | 22.42 | 21.61 | 24.08 | 21.35 |
| `benchMany_cxx17::6 funcs` | 11.12 | 10.80 | 8.03 | 8.06 |
| `benchMany_cxx17::6 funcs (non-many)` | 27.63 | 38.20 | 29.35 | 25.40 |
| `benchMany_cxx17::7 funcs` | 11.72 | 11.98 | 8.37 | 8.17 |
| `benchMany_cxx17::7 funcs (non-many)` | 32.16 | 32.87 | 33.97 | 30.55 |
| `benchMany_cxx17::8 funcs` | 10.95 | 11.20 | 8.19 | 8.02 |
| `benchMany_cxx17::8 funcs (non-many)` | 35.75 | 35.02 | 38.54 | 35.34 |
| `benchMany_cxx17::9 funcs` | 14.43 | 14.23 | 11.13 | 12.79 |
| `benchMany_cxx17::9 funcs (non-many)` | 56.72 | 39.07 | 43.00 | 39.15 |
| `benchMany_cxx17::10 funcs` | 14.48 | 14.71 | 11.21 | 13.09 |
| `benchMany_cxx17::10 funcs (non-many)` | 61.81 | 43.63 | 47.96 | 43.26 |
| `benchMany_cxx17::11 funcs` | 15.98 | 16.20 | 10.92 | 11.15 |
| `benchMany_cxx17::11 funcs (non-many)` | 49.04 | 62.47 | 52.29 | 46.57 |
| `benchMany_cxx17::12 funcs` | 14.90 | 14.55 | 11.31 | 13.30 |
| `benchMany_cxx17::12 funcs (non-many)` | 75.61 | 52.43 | 62.01 | 52.54 |
| `benchMany_cxx17::13 funcs` | 17.18 | 16.88 | 13.75 | 14.38 |
| `benchMany_cxx17::13 funcs (non-many)` | 79.11 | 73.92 | 66.97 | 56.90 |
| `benchMany_cxx17::14 funcs` | 18.24 | 17.68 | 13.91 | 17.04 |
| `benchMany_cxx17::14 funcs (non-many)` | 68.22 | 65.28 | 72.19 | 61.38 |
| `benchMany_cxx17::15 funcs` | 19.08 | 19.08 | 13.39 | 15.29 |
| `benchMany_cxx17::15 funcs (non-many)` | 72.41 | 88.76 | 82.62 | 65.63 |
| `benchMany_cxx17::16 funcs` | 17.25 | 16.89 | 13.41 | 16.59 |
| `benchMany_cxx17::16 funcs (non-many)` | 102.22 | 73.21 | 88.71 | 69.96 |
| `benchMany_cxx17::8 funcs (mixed domains)` | 11.11 | 11.18 | 8.92 | 8.35 |
| `benchMany_cxx17::8 funcs (mixed domains, non-many)` | 35.78 | 35.00 | 38.51 | 33.82 |
| `benchMany_cxx17::16 funcs (mixed domains)` | 17.21 | 16.90 | 13.75 | 15.82 |
| `benchMany_cxx17::16 funcs (mixed domains, non-many)` | 102.29 | 73.09 | 88.19 | 70.00 |
| `benchMany_cxx17::8 funcs (complex outputs, bulk)` | 38.53 | 38.53 | 43.56 | 39.48 |
| `benchMany_cxx20::1 funcs` | 7.76 | 7.54 | 5.59 | 5.25 |
| `benchMany_cxx20::1 funcs (non-many)` | 4.58 | 4.39 | 5.24 | 4.62 |
| `benchMany_cxx20::2 funcs` | 8.53 | 7.85 | 5.65 | 5.26 |
| `benchMany_cxx20::2 funcs (non-many)` | 9.38 | 9.28 | 10.13 | 8.55 |
| `benchMany_cxx20::3 funcs` | 8.52 | 8.65 | 6.14 | 5.08 |
| `benchMany_cxx20::3 funcs (non-many)` | 19.27 | 13.42 | 15.11 | 12.87 |
| `benchMany_cxx20::4 funcs` | 8.27 | 8.99 | 5.60 | 5.34 |
| `benchMany_cxx20::4 funcs (non-many)` | 18.43 | 17.98 | 19.62 | 17.41 |
| `benchMany_cxx20::5 funcs` | 14.21 | 10.91 | 8.13 | 7.98 |
| `benchMany_cxx20::5 funcs (non-many)` | 22.43 | 21.52 | 24.33 | 21.21 |
| `benchMany_cxx20::6 funcs` | 14.35 | 10.95 | 8.33 | 8.09 |
| `benchMany_cxx20::6 funcs (non-many)` | 27.66 | 38.19 | 29.35 | 25.63 |
| `benchMany_cxx20::7 funcs` | 11.73 | 11.79 | 8.33 | 8.60 |
| `benchMany_cxx20::7 funcs (non-many)` | 32.17 | 32.79 | 33.35 | 29.73 |
| `benchMany_cxx20::8 funcs` | 10.95 | 11.01 | 8.41 | 8.38 |
| `benchMany_cxx20::8 funcs (non-many)` | 35.73 | 35.03 | 38.62 | 33.78 |
| `benchMany_cxx20::9 funcs` | 14.32 | 14.34 | 11.16 | 12.84 |
| `benchMany_cxx20::9 funcs (non-many)` | 56.76 | 39.08 | 43.32 | 38.08 |
| `benchMany_cxx20::10 funcs` | 14.90 | 14.75 | 11.44 | 13.11 |
| `benchMany_cxx20::10 funcs (non-many)` | 61.77 | 43.66 | 48.09 | 42.17 |
| `benchMany_cxx20::11 funcs` | 16.15 | 16.21 | 11.63 | 11.19 |
| `benchMany_cxx20::11 funcs (non-many)` | 49.00 | 62.52 | 52.59 | 48.55 |
| `benchMany_cxx20::12 funcs` | 15.08 | 15.29 | 9.78 | 12.77 |
| `benchMany_cxx20::12 funcs (non-many)` | 75.59 | 52.44 | 61.87 | 55.03 |
| `benchMany_cxx20::13 funcs` | 17.03 | 16.97 | 13.84 | 14.20 |
| `benchMany_cxx20::13 funcs (non-many)` | 79.10 | 73.91 | 66.98 | 56.88 |
| `benchMany_cxx20::14 funcs` | 17.71 | 17.61 | 13.76 | 17.01 |
| `benchMany_cxx20::14 funcs (non-many)` | 68.26 | 65.94 | 73.18 | 64.94 |
| `benchMany_cxx20::15 funcs` | 18.92 | 18.85 | 13.60 | 15.36 |
| `benchMany_cxx20::15 funcs (non-many)` | 72.24 | 88.66 | 81.99 | 67.72 |
| `benchMany_cxx20::16 funcs` | 17.45 | 16.89 | 14.58 | 16.33 |
| `benchMany_cxx20::16 funcs (non-many)` | 102.19 | 73.21 | 88.36 | 72.09 |
| `benchMany_cxx20::8 funcs (mixed domains)` | 10.93 | 10.76 | 8.14 | 7.96 |
| `benchMany_cxx20::8 funcs (mixed domains, non-many)` | 35.74 | 34.99 | 38.75 | 34.06 |
| `benchMany_cxx20::16 funcs (mixed domains)` | 17.24 | 16.89 | 13.82 | 15.91 |
| `benchMany_cxx20::16 funcs (mixed domains, non-many)` | 102.13 | 73.10 | 87.70 | 72.14 |
| `benchMany_cxx20::8 funcs (complex outputs, bulk)` | 38.47 | 38.54 | 43.33 | 39.03 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 build` | 61165.15 | 60677.76 | 78371.70 | 82089.45 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16 build` | 88934.38 | 88077.66 | 106211.95 | 107031.89 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16 build` | 1522891.43 | 1504598.25 | 1951045.83 | 2034865.83 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8 build` | 239395.25 | 236470.05 | 347781.45 | 389170.55 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8 build` | 315753.51 | 312675.85 | 474028.50 | 523252.45 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8 build` | 2540961.60 | 2568958.50 | 3827850.00 | 4246491.33 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8 build` | 3336448.50 | 3396934.00 | 5137245.00 | 5455129.00 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16` | 51.63 | 52.19 | 900.62 | 1303.41 |
| `bench_ND_cxx17::F:ℝ²→ℝ³, D=16` | 112.16 | 112.31 | 998.02 | 1157.60 |
| `bench_ND_cxx17::F:ℝ³→ℝ², D=16` | 1272.54 | 1279.32 | 14840.77 | 18485.69 |
| `bench_ND_cxx17::F:ℝ³→ℝ³, D=8` | 197.04 | 197.25 | 2104.28 | 2788.58 |
| `bench_ND_cxx17::F:ℝ³→ℝ⁴, D=8` | 107.54 | 107.10 | 1960.59 | 2142.58 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ³, D=8` | 2606.04 | 2606.00 | 17271.14 | 20893.69 |
| `bench_ND_cxx17::F:ℝ⁴→ℝ⁴, D=8` | 1725.43 | 1718.02 | 16136.95 | 19258.33 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic point` | 51.58 | 52.24 | 900.00 | 1298.85 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical point` | 51.58 | 52.21 | 899.64 | 1299.57 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 array` | 51.59 | 52.22 | 901.71 | 1304.36 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 variadic` | 92.96 | 92.82 | 902.45 | 1302.82 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 manual loop` | 52823.32 | 52897.87 | 922685.68 | 1336022.94 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 canonical batch` | 113119.75 | 111943.04 | 922154.78 | 1307962.56 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 generic loop` | 52782.89 | 52845.80 | 922462.88 | 1340326.31 |
| `bench_ND_cxx17::F:ℝ²→ℝ², D=16 container batch` | 112508.64 | 113545.81 | 920775.26 | 1328822.12 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 build` | 60955.09 | 60889.15 | 78682.84 | 82974.69 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16 build` | 88646.23 | 88674.36 | 106259.39 | 106704.81 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16 build` | 1523178.50 | 1529561.86 | 1972320.33 | 2054643.80 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8 build` | 232991.29 | 233790.06 | 346429.54 | 388536.15 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8 build` | 307652.46 | 309811.66 | 473820.92 | 522542.30 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8 build` | 2483787.00 | 2535113.50 | 3829209.33 | 4037419.67 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8 build` | 3265262.33 | 3366031.67 | 5149277.50 | 5661660.50 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16` | 51.58 | 52.28 | 901.56 | 1298.55 |
| `bench_ND_cxx20::F:ℝ²→ℝ³, D=16` | 111.93 | 112.32 | 1003.14 | 1179.47 |
| `bench_ND_cxx20::F:ℝ³→ℝ², D=16` | 1263.41 | 1272.84 | 14834.10 | 18165.06 |
| `bench_ND_cxx20::F:ℝ³→ℝ³, D=8` | 197.54 | 197.48 | 2141.06 | 2620.71 |
| `bench_ND_cxx20::F:ℝ³→ℝ⁴, D=8` | 104.69 | 104.74 | 1983.74 | 2168.79 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ³, D=8` | 2603.89 | 2606.64 | 17253.59 | 21276.61 |
| `bench_ND_cxx20::F:ℝ⁴→ℝ⁴, D=8` | 1727.75 | 1718.43 | 16116.72 | 19273.11 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic point` | 51.61 | 52.22 | 901.35 | 1309.14 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical point` | 51.56 | 52.21 | 902.43 | 1302.30 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 array` | 51.58 | 52.20 | 901.63 | 1300.36 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 variadic` | 93.17 | 97.05 | 901.91 | 1299.52 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 manual loop` | 52875.82 | 52896.79 | 923805.52 | 1343363.56 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 canonical batch` | 111150.31 | 111804.50 | 925206.13 | 1304838.83 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 span batch` | 92218.00 | 92261.09 | 923622.74 | 1307846.17 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 generic loop` | 52855.08 | 52830.86 | 925420.32 | 1339538.44 |
| `bench_ND_cxx20::F:ℝ²→ℝ², D=16 container batch` | 112525.44 | 113781.34 | 924652.28 | 1337858.24 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.16 | 2.19 | 3.40 | 3.53 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.35 | 0.35 | 0.35 | 0.34 |
| `bench_horner_cxx17::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.35 | 0.36 | 0.35 | 0.34 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.37 | 4.38 | 5.15 | 4.67 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.66 | 0.67 | 0.67 | 0.71 |
| `bench_horner_cxx17::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.66 | 0.67 | 0.68 | 0.71 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.56 | 7.54 | 7.63 | 6.80 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 1.05 | 1.05 | 1.05 | 1.14 |
| `bench_horner_cxx17::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 1.05 | 1.05 | 1.05 | 1.13 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 11.12 | 11.18 | 10.56 | 10.21 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.43 | 1.45 | 1.43 | 1.56 |
| `bench_horner_cxx17::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.43 | 1.44 | 1.40 | 1.56 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.53 | 2.45 | 5.48 | 5.56 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.81 | 2.54 | 4.60 | 4.28 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.57 | 2.29 | 4.02 | 3.88 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.37 | 2.08 | 3.86 | 3.65 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 4.76 | 4.65 | 8.43 | 7.34 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 5.18 | 4.66 | 7.61 | 5.79 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.85 | 4.49 | 6.49 | 5.16 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.48 | 4.14 | 6.16 | 4.94 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 7.61 | 8.01 | 11.23 | 10.30 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 8.24 | 7.67 | 10.67 | 8.37 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 7.98 | 7.56 | 9.60 | 7.67 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 7.49 | 7.12 | 9.09 | 7.41 |
| `bench_horner_cxx17::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 12.57 | 12.71 | 15.20 | 12.78 |
| `bench_horner_cxx17::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 12.25 | 11.85 | 13.76 | 11.92 |
| `bench_horner_cxx17::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 11.79 | 11.52 | 12.76 | 10.82 |
| `bench_horner_cxx17::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 11.25 | 10.57 | 12.21 | 10.90 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 8.86 | 8.79 | 6.24 | 6.36 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.86 | 3.11 | 2.46 | 2.93 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 2.47 | 2.36 | 2.89 | 3.23 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 2.28 | 2.03 | 2.06 | 2.21 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 18.66 | 18.66 | 13.86 | 13.38 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 6.71 | 6.66 | 6.53 | 7.08 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 4.46 | 4.27 | 5.38 | 6.01 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 4.38 | 4.31 | 4.40 | 4.58 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 28.07 | 27.99 | 21.82 | 21.30 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 10.32 | 9.91 | 8.91 | 11.15 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 8.59 | 8.27 | 9.17 | 10.06 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 7.04 | 6.74 | 6.97 | 7.45 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 37.41 | 37.32 | 30.03 | 30.68 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 14.61 | 14.30 | 13.26 | 13.55 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 12.33 | 11.82 | 12.11 | 13.43 |
| `bench_horner_cxx17::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 9.37 | 9.09 | 9.54 | 10.47 |
| `bench_horner_cxx17::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 8.57 | 8.09 | 7.05 | 7.23 |
| `bench_horner_cxx17::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 27.55 | 27.41 | 16.90 | 16.90 |
| `bench_horner_cxx17::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 132.68 | 132.33 | 63.14 | 61.37 |
| `bench_horner_cxx17::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 337.75 | 332.41 | 163.73 | 167.71 |
| `bench_horner_cxx17::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 35.65 | 35.60 | 35.48 | 39.77 |
| `bench_horner_cxx17::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 236.25 | 232.58 | 202.93 | 238.80 |
| `bench_horner_cxx17::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2171.48 | 2115.59 | 1764.77 | 1818.28 |
| `bench_horner_cxx17::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 146.50 | 136.81 | 129.25 | 128.10 |
| `bench_horner_cxx17::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 2003.89 | 2188.22 | 1592.92 | 1614.15 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=No, scalar-runtime` | 2.15 | 2.17 | 3.37 | 3.47 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-aligned` | 0.35 | 0.35 | 0.34 | 0.34 |
| `bench_horner_cxx20::Dim=1, nCoeffs=8, SIMD=Yes, simd-unaligned` | 0.35 | 0.36 | 0.35 | 0.34 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=No, scalar-runtime` | 4.36 | 4.39 | 5.15 | 4.65 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-aligned` | 0.67 | 0.66 | 0.68 | 0.71 |
| `bench_horner_cxx20::Dim=1, nCoeffs=16, SIMD=Yes, simd-unaligned` | 0.66 | 0.67 | 0.67 | 0.71 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=No, scalar-runtime` | 7.56 | 7.56 | 7.63 | 6.75 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-aligned` | 1.06 | 1.05 | 1.06 | 1.13 |
| `bench_horner_cxx20::Dim=1, nCoeffs=24, SIMD=Yes, simd-unaligned` | 1.05 | 1.05 | 1.06 | 1.13 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=No, scalar-runtime` | 11.11 | 11.18 | 10.56 | 10.21 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-aligned` | 1.44 | 1.44 | 1.40 | 1.56 |
| `bench_horner_cxx20::Dim=1, nCoeffs=32, SIMD=Yes, simd-unaligned` | 1.43 | 1.44 | 1.43 | 1.56 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=8, SIMD=No, many-runtime` | 2.77 | 2.46 | 5.36 | 5.52 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=8, SIMD=No, many-runtime` | 2.79 | 2.53 | 4.64 | 4.29 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=8, SIMD=No, many-runtime` | 2.59 | 2.33 | 4.10 | 3.90 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=8, SIMD=No, many-runtime` | 2.31 | 2.07 | 3.90 | 3.67 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=16, SIMD=No, many-runtime` | 4.88 | 4.65 | 8.47 | 7.32 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=16, SIMD=No, many-runtime` | 5.14 | 4.67 | 7.66 | 5.80 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=16, SIMD=No, many-runtime` | 4.96 | 4.51 | 6.45 | 5.19 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=16, SIMD=No, many-runtime` | 4.47 | 4.14 | 6.09 | 4.86 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=24, SIMD=No, many-runtime` | 8.14 | 8.01 | 11.63 | 10.30 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=24, SIMD=No, many-runtime` | 8.28 | 7.69 | 10.67 | 8.46 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=24, SIMD=No, many-runtime` | 7.97 | 7.51 | 9.61 | 7.66 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=24, SIMD=No, many-runtime` | 7.47 | 7.12 | 9.09 | 7.72 |
| `bench_horner_cxx20::Dim=1, horner_many M=4, nCoeffs=32, SIMD=No, many-runtime` | 12.76 | 12.72 | 15.03 | 12.73 |
| `bench_horner_cxx20::Dim=1, horner_many M=8, nCoeffs=32, SIMD=No, many-runtime` | 12.19 | 11.88 | 13.68 | 11.84 |
| `bench_horner_cxx20::Dim=1, horner_many M=12, nCoeffs=32, SIMD=No, many-runtime` | 11.87 | 11.51 | 12.78 | 10.72 |
| `bench_horner_cxx20::Dim=1, horner_many M=16, nCoeffs=32, SIMD=No, many-runtime` | 11.34 | 10.68 | 12.40 | 10.98 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=8` | 8.86 | 8.79 | 6.12 | 6.36 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=8` | 2.91 | 3.10 | 2.61 | 2.91 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=8` | 2.46 | 2.35 | 2.87 | 3.26 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=8` | 2.20 | 2.03 | 2.03 | 2.20 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=16` | 18.66 | 18.66 | 13.12 | 13.40 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=16` | 6.66 | 6.62 | 6.44 | 7.07 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=16` | 4.49 | 4.22 | 5.36 | 6.02 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=16` | 4.39 | 4.31 | 4.44 | 4.57 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=24` | 28.01 | 27.98 | 20.63 | 21.65 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=24` | 10.31 | 9.91 | 8.96 | 11.14 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=24` | 8.60 | 8.27 | 9.15 | 10.04 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=24` | 6.96 | 6.76 | 6.92 | 7.45 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=4, nCoeffs=32` | 37.32 | 37.36 | 28.52 | 30.67 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=8, nCoeffs=32` | 14.61 | 14.26 | 13.27 | 13.55 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=12, nCoeffs=32` | 12.34 | 11.82 | 12.21 | 13.41 |
| `bench_horner_cxx20::Dim=1, horner_transposed_scalar M=16, nCoeffs=32` | 9.37 | 9.09 | 9.57 | 10.47 |
| `bench_horner_cxx20::Dim=2, nCoeffs=4, SIMD=No, ND2D` | 8.64 | 8.09 | 7.05 | 7.39 |
| `bench_horner_cxx20::Dim=2, nCoeffs=8, SIMD=No, ND2D` | 27.55 | 27.36 | 16.82 | 19.00 |
| `bench_horner_cxx20::Dim=2, nCoeffs=16, SIMD=No, ND2D` | 132.67 | 132.23 | 64.50 | 64.04 |
| `bench_horner_cxx20::Dim=2, nCoeffs=24, SIMD=No, ND2D` | 337.68 | 332.38 | 172.76 | 167.92 |
| `bench_horner_cxx20::Dim=3, nCoeffs=4, SIMD=No, ND3D` | 35.50 | 35.60 | 35.48 | 40.11 |
| `bench_horner_cxx20::Dim=3, nCoeffs=8, SIMD=No, ND3D` | 236.22 | 232.66 | 201.19 | 237.79 |
| `bench_horner_cxx20::Dim=3, nCoeffs=16, SIMD=No, ND3D` | 2104.33 | 2135.20 | 1769.70 | 1821.11 |
| `bench_horner_cxx20::Dim=4, nCoeffs=4, SIMD=No, ND4D` | 146.13 | 137.09 | 126.20 | 129.89 |
| `bench_horner_cxx20::Dim=4, nCoeffs=8, SIMD=No, ND4D` | 2148.63 | 2074.71 | 1584.29 | 1559.02 |

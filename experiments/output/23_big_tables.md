# Large tables: 2 chains, 500 warm-up + 500 draws, limit 900 s (largest-column margin coordinates)

Every fit evaluates the prior with a dense RC x (RC-1) basis; that product dominates the time per gradient from about 24x24.

| size   |   cells | param            | status                                 |   seconds |   leapfrogs |   ms_per_gradient |   min_ess |   max_rhat |   ess_per_s |
|:-------|--------:|:-----------------|:---------------------------------------|----------:|------------:|------------------:|----------:|-----------:|------------:|
| 10x10  |     100 | odds-ratio logit | ok                                     |      1.25 |       15    |              0.09 |   1377.99 |       1.03 |     1105.87 |
| 10x10  |     100 | position logit   | ok                                     |      2.06 |       20.6  |              0.09 |    731    |       1.01 |      354.58 |
| 10x10  |     100 | adjusted row     | ok                                     |      2.38 |       12.13 |              0.15 |    201.72 |       1.01 |       84.87 |
| 10x10  |     100 | ILR + log volume | ok                                     |      4.97 |       62.97 |              0.04 |    904.46 |       1.01 |      181.96 |
| 10x10  |     100 | log cells        | ok                                     |      3.47 |      118.42 |              0.03 |   1186.63 |       1.02 |      341.96 |
| 10x10  |     100 | adjusted table   | ok                                     |     51.63 |       52.81 |              1.15 |    378.43 |       1.01 |        7.33 |
| 5x40   |     200 | odds-ratio logit | ok                                     |      3.13 |       15    |              0.21 |   1242.1  |       1.02 |      397.44 |
| 5x40   |     200 | position logit   | ok                                     |     11.11 |       47    |              0.2  |    755.23 |       1.02 |       67.99 |
| 5x40   |     200 | adjusted row     | ok                                     |      4.6  |       18.92 |              0.26 |     74.61 |       1.04 |       16.23 |
| 5x40   |     200 | ILR + log volume | ok                                     |      9.81 |       31    |              0.12 |    777.55 |       1.02 |       79.3  |
| 5x40   |     200 | log cells        | ok                                     |      5.77 |       63    |              0.09 |   1350.88 |       1.02 |      233.94 |
| 5x40   |     200 | adjusted table   | ok                                     |    384.41 |       80.98 |              5.12 |    319.32 |       1.02 |        0.83 |
| 16x16  |     256 | odds-ratio logit | ok                                     |      5.46 |       15    |              0.36 |   1247.85 |       1.02 |      228.67 |
| 16x16  |     256 | position logit   | ok                                     |      9.73 |       30.87 |              0.32 |    926.07 |       1.01 |       95.21 |
| 16x16  |     256 | adjusted row     | ok                                     |      8.97 |       15.29 |              0.4  |     73.78 |       1.03 |        8.22 |
| 16x16  |     256 | ILR + log volume | ok                                     |     27.57 |       62.87 |              0.22 |   1172.98 |       1.02 |       42.54 |
| 16x16  |     256 | log cells        | ok                                     |     17.15 |       75.22 |              0.2  |    811.17 |       1.02 |       47.31 |
| 16x16  |     256 | adjusted table   | TimeoutError                           |    900.08 |      nan    |            nan    |    nan    |     nan    |      nan    |
| 24x24  |     576 | odds-ratio logit | ok                                     |     20.85 |       15    |              1.26 |    843.1  |       1.02 |       40.44 |
| 24x24  |     576 | position logit   | ok                                     |     70.97 |       63    |              1.17 |    948.17 |       1.02 |       13.36 |
| 24x24  |     576 | adjusted row     | ok                                     |     32.31 |       15    |              1.4  |     54.04 |       1.03 |        1.67 |
| 24x24  |     576 | ILR + log volume | ok                                     |    139.9  |       62.9  |              1.09 |   1187.64 |       1.02 |        8.49 |
| 24x24  |     576 | log cells        | ok                                     |     64.4  |       63    |              0.87 |    951.94 |       1.02 |       14.78 |
| 24x24  |     576 | adjusted table   | TimeoutError                           |    900.36 |      nan    |            nan    |    nan    |     nan    |      nan    |
| 32x32  |    1024 | odds-ratio logit | ok                                     |    192.07 |       15    |             10.74 |    766.07 |       1.01 |        3.99 |
| 32x32  |    1024 | position logit   | ok                                     |    505.03 |       47    |             10.29 |    585.88 |       1.01 |        1.16 |
| 32x32  |    1024 | adjusted row     | not converged                          |    237.05 |       15    |             10.38 |     13.84 |       1.13 |        0.06 |
| 32x32  |    1024 | ILR + log volume | TimeoutError                           |    900.93 |      nan    |            nan    |    nan    |     nan    |      nan    |
| 32x32  |    1024 | log cells        | ok                                     |    593.4  |       63    |              8.37 |   1028.28 |       1.03 |        1.73 |
| 32x32  |    1024 | adjusted table   | not run (timed out at 16x16 and 24x24) |    nan    |      nan    |            nan    |    nan    |     nan    |      nan    |

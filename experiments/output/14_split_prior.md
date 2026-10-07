# Prior on interactions alone against equal sds, split basis, Scottish grand total (eps 0.1)

|                                                          |   failed |   max_rhat |   min_ess |   leapfrogs |   ess_per_1000_grad |   div |   seconds |
|:---------------------------------------------------------|---------:|-----------:|----------:|------------:|--------------------:|------:|----------:|
| ('indep', 1.0, 'position logit', 'equal sd')             |        0 |       1.01 |    560.55 |       15.86 |               35.34 |     0 |      0.43 |
| ('indep', 1.0, 'odds-ratio logit', 'equal sd')           |        0 |       1    |    979.42 |       15.11 |               64.83 |     0 |      0.4  |
| ('indep', 1.0, 'adjusted row', 'equal sd')               |        0 |       1.01 |    491.15 |       22.16 |               22.02 |     0 |      0.79 |
| ('indep', 1.0, 'adjusted table', 'equal sd')             |        0 |       1.01 |    325.31 |       38.03 |                8.55 |     0 |      5.46 |
| ('indep', 1.0, 'position logit', 'interactions only')    |        0 |       1.01 |    528.95 |       15.67 |               33.93 |     0 |      0.44 |
| ('indep', 1.0, 'odds-ratio logit', 'interactions only')  |        0 |       1    |    584.65 |       15.39 |               37.98 |     0 |      0.37 |
| ('indep', 1.0, 'adjusted row', 'interactions only')      |        0 |       1    |    625.07 |       26.28 |               23.82 |     0 |      1    |
| ('indep', 1.0, 'adjusted table', 'interactions only')    |        0 |       1.01 |    409.28 |       32.88 |               12.43 |     0 |      5.64 |
| ('indep', 2.0, 'position logit', 'equal sd')             |        0 |       1.01 |    448.22 |      145.1  |                3.26 |     0 |      3.06 |
| ('indep', 2.0, 'odds-ratio logit', 'equal sd')           |        0 |       1.02 |    211.45 |       33.12 |                6.45 |     2 |      0.99 |
| ('indep', 2.0, 'adjusted row', 'equal sd')               |        0 |       1.01 |    322.22 |       22.73 |               14.07 |     0 |      1.17 |
| ('indep', 2.0, 'adjusted table', 'equal sd')             |        0 |       1.03 |    144.84 |       32.12 |                4.54 |     0 |      7.07 |
| ('indep', 2.0, 'position logit', 'interactions only')    |        0 |       1.02 |    234.39 |      215.14 |                1.11 |     6 |      5.05 |
| ('indep', 2.0, 'odds-ratio logit', 'interactions only')  |        0 |       1.02 |    268.84 |       48.13 |                5.61 |     3 |      1.62 |
| ('indep', 2.0, 'adjusted row', 'interactions only')      |        0 |       1.01 |    303.79 |       23.16 |               13.11 |     0 |      1.24 |
| ('indep', 2.0, 'adjusted table', 'interactions only')    |        0 |       1.01 |    406.43 |       38.79 |               10.38 |     0 |      6.68 |
| ('actual', 1.0, 'position logit', 'equal sd')            |        0 |       1.01 |    593.73 |       56.97 |               10.42 |     0 |      1.18 |
| ('actual', 1.0, 'odds-ratio logit', 'equal sd')          |        0 |       1.01 |    456.3  |       29.51 |               15.6  |     0 |      0.61 |
| ('actual', 1.0, 'adjusted row', 'equal sd')              |        0 |       1.02 |    383.87 |       19.68 |               19.52 |     0 |      0.81 |
| ('actual', 1.0, 'adjusted table', 'equal sd')            |        0 |       1.03 |    200.5  |       31.39 |                6.39 |     0 |      5.97 |
| ('actual', 1.0, 'position logit', 'interactions only')   |        0 |       1.02 |    200.56 |       69.27 |                2.91 |     0 |      1.97 |
| ('actual', 1.0, 'odds-ratio logit', 'interactions only') |        0 |       1.01 |    301.1  |       29.04 |               10.35 |     0 |      0.76 |
| ('actual', 1.0, 'adjusted row', 'interactions only')     |        0 |       1.01 |    616.92 |       27.49 |               22.43 |     0 |      0.97 |
| ('actual', 1.0, 'adjusted table', 'interactions only')   |        0 |       1.02 |    209.36 |       25.54 |                8.05 |     0 |      5.01 |
| ('actual', 2.0, 'position logit', 'equal sd')            |        0 |       1.01 |    237.62 |      772.26 |                0.31 |     1 |     13.78 |
| ('actual', 2.0, 'odds-ratio logit', 'equal sd')          |        0 |       1.02 |    235.35 |       98.02 |                2.38 |     1 |      2.55 |
| ('actual', 2.0, 'adjusted row', 'equal sd')              |        0 |       1.01 |    230.57 |       19.41 |               11.86 |     0 |      1.02 |
| ('actual', 2.0, 'adjusted table', 'equal sd')            |        0 |       1.02 |    116.01 |       36.18 |                3.21 |     0 |      7.72 |
| ('actual', 2.0, 'position logit', 'interactions only')   |        1 |       1.07 |     67.43 |     1011.94 |                0.07 |     7 |     18.05 |
| ('actual', 2.0, 'odds-ratio logit', 'interactions only') |        0 |       1.02 |    218.81 |      223.38 |                0.97 |     0 |      5.2  |
| ('actual', 2.0, 'adjusted row', 'interactions only')     |        0 |       1.01 |    201.56 |       24.26 |                8.3  |     0 |      1.38 |
| ('actual', 2.0, 'adjusted table', 'interactions only')   |        0 |       1.01 |    388.69 |       40.79 |                9.55 |     0 |      6.92 |

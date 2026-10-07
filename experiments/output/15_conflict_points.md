# Conflict points: median ESS per 1000 gradients over three areas (converged runs)

|                                        |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:---------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| (0.5, 'indep', 'equal sd')             |                6.5 |         6.1 |             86.9 |              129.7 |           42.9 |             29.1 |
| (0.5, 'raked', 'equal sd')             |                3.2 |         8.7 |             60.5 |               84.4 |           43.9 |             18   |
| (0.5, 'national', 'equal sd')          |                2.7 |         4.8 |             72.4 |               58.6 |           37.1 |             16.4 |
| (0.5, 'national', 'interactions only') |                2   |         5.5 |             70.5 |               59.4 |           40.9 |             30   |
| (2.0, 'indep', 'equal sd')             |                0.5 |         0.4 |              5.4 |               13.8 |           25.6 |             13.6 |
| (2.0, 'raked', 'equal sd')             |                0.4 |         0.2 |              2.8 |                6.5 |            8.7 |             10.6 |
| (2.0, 'national', 'equal sd')          |                0.3 |         0.4 |              2   |                3.4 |           13.2 |             12.9 |
| (2.0, 'national', 'interactions only') |                0.3 |         0.1 |              0.9 |                2.6 |           15.4 |             19.6 |

## Part B: volume (median area, raked centre, sigma_b 1), ESS per 1000 gradients

|              N |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|---------------:|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
|   296.92       |              17.06 |       17.83 |            21.74 |              34.88 |          26.04 |            20.14 |
| 29692          |               1.21 |        0.83 |            44    |              46.4  |          25.47 |            17.42 |
|     2.9692e+06 |               0.01 |        0.13 |            44.04 |              29.54 |          34.52 |            14.64 |

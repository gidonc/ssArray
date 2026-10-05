# The single-area model in six coordinate systems

## Part 1: the stated Jacobian against finite differences

Six random points per row. 'Offset' is the mean of (finite difference - stated); it should be zero, or a constant for the ILR parameterisation, whose Jacobian is stated up to a constant. 'Spread' is the largest deviation from that mean and should be near zero.

| Size | Parameterisation | Offset | Spread |
|---|---|---|---|
| 3x3 | ILR + log volume | 1.10e+00 | 9.3e-11 |
| 3x3 | log cells | -1.87e-10 | 1.7e-11 |
| 3x3 | between bounds | 2.06e-09 | 9.3e-09 |
| 3x3 | cell logit | 1.63e-09 | 5.8e-09 |
| 3x3 | adjusted row | -9.34e-10 | 6.4e-09 |
| 3x3 | adjusted table | 6.35e-09 | 2.9e-08 |
| 7x7 | ILR + log volume | 1.95e+00 | 2.6e-10 |
| 7x7 | log cells | -3.75e-10 | 9.2e-11 |
| 7x7 | between bounds | -3.33e-09 | 7.8e-09 |
| 7x7 | cell logit | 9.01e-10 | 5.2e-09 |
| 7x7 | adjusted row | 1.91e-09 | 9.5e-09 |
| 7x7 | adjusted table | 1.58e-08 | 3.1e-08 |

## Part 2: sampling the same posterior

2 chains, 500 warm-up + 500 draws each, started at the independence table. Placeholder prior: ILR coordinates ~ normal(0, 2). Summaries are of the log cells. 'Gap' is the largest difference in a posterior mean from the adjusted-row run, in posterior sds; 'sd ratio' is the range of the ratio of posterior sds to that run. ESS is the smallest over cells, of 1,000 draws.

| Size | eps | Parameterisation | Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | Seconds | Gap | sd ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3x3 | 1 | ILR + log volume | 0.00845 | 420 | 0 | 1% | 144 | 1.007 | 0.34 | 1 | 0.08 | 0.93-1.07 |
| 3x3 | 1 | log cells | 0.00727 | 457 | 0 | 5% | 83 | 1.027 | 0.18 | 1 | 0.19 | 0.87-1.18 |
| 3x3 | 1 | between bounds | 0.186 | 20 | 0 | 0% | 326 | 1.008 | 16.27 | 0 | 0.08 | 1.00-1.09 |
| 3x3 | 1 | cell logit | 0.302 | 13 | 0 | 0% | 386 | 1.006 | 29.29 | 0 | 0.16 | 0.84-1.23 |
| 3x3 | 1 | adjusted row | 0.373 | 11 | 0 | 0% | 387 | 1.006 | 34.23 | 0 | 0.00 | 1.00-1.00 |
| 3x3 | 1 | adjusted table | 0.349 | 12 | 0 | 0% | 575 | 1.005 | 47.00 | 0 | 0.06 | 0.88-1.06 |
| 3x3 | 0.1 | ILR + log volume | 0.000777 | 894 | 0 | 82% | 9 | 1.478 | 0.01 | 3 | 0.40 | 0.55-1.16 |
| 3x3 | 0.1 | log cells | 0.00163 | 790 | 83 | 67% | 3 | 1.254 | 0.00 | 2 | 1.17 | 0.22-1.20 |
| 3x3 | 0.1 | between bounds | 0.263 | 15 | 0 | 0% | 475 | 1.007 | 32.73 | 0 | 0.07 | 0.90-1.12 |
| 3x3 | 0.1 | cell logit | 0.275 | 14 | 0 | 0% | 286 | 1.008 | 20.28 | 0 | 0.08 | 0.92-1.20 |
| 3x3 | 0.1 | adjusted row | 0.367 | 12 | 0 | 0% | 514 | 1.001 | 44.17 | 0 | 0.00 | 1.00-1.00 |
| 3x3 | 0.1 | adjusted table | 0.34 | 12 | 0 | 0% | 586 | 1.003 | 47.21 | 0 | 0.08 | 0.97-1.02 |
| 7x7 | 1 | ILR + log volume | 0.0109 | 485 | 8 | 0% | 116 | 1.019 | 0.24 | 6 | 0.17 | 0.25-1.07 |
| 7x7 | 1 | log cells | 0.0112 | 474 | 0 | 0% | 61 | 1.036 | 0.13 | 4 | 0.18 | 0.17-1.10 |
| 7x7 | 1 | between bounds | 0.0101 | 390 | 2 | 0% | 48 | 1.044 | 0.12 | 11 | 0.43 | 0.90-1.69 |
| 7x7 | 1 | cell logit | 0.0661 | 66 | 0 | 0% | 34 | 1.052 | 0.50 | 2 | 0.33 | 0.91-1.63 |
| 7x7 | 1 | adjusted row | 0.17 | 31 | 0 | 0% | 54 | 1.030 | 1.76 | 2 | 0.00 | 1.00-1.00 |
| 7x7 | 1 | adjusted table | 0.159 | 31 | 0 | 0% | 146 | 1.020 | 4.75 | 7 | 0.48 | 0.93-1.79 |
| 7x7 | 0.1 | ILR + log volume | 0.00124 | 1023 | 0 | 100% | 20 | 1.131 | 0.02 | 11 | 0.53 | 0.18-1.14 |
| 7x7 | 0.1 | log cells | 0.00129 | 1023 | 0 | 100% | 17 | 1.082 | 0.02 | 8 | 0.40 | 0.05-1.19 |
| 7x7 | 0.1 | between bounds | 0.0155 | 346 | 1 | 0% | 47 | 1.040 | 0.14 | 9 | 0.14 | 0.81-1.15 |
| 7x7 | 0.1 | cell logit | 0.0457 | 93 | 0 | 0% | 30 | 1.062 | 0.32 | 2 | 0.10 | 0.81-1.18 |
| 7x7 | 0.1 | adjusted row | 0.171 | 30 | 0 | 0% | 65 | 1.037 | 2.15 | 2 | 0.00 | 1.00-1.00 |
| 7x7 | 0.1 | adjusted table | 0.158 | 31 | 0 | 0% | 181 | 1.008 | 5.85 | 8 | 0.13 | 0.91-1.23 |

## Part 3: longer runs of the sequential parameterisations at 7x7

2 chains, 1,000 warm-up + 3,000 draws each, eps = 1. Same columns; ESS is of 6,000 draws; gap and sd ratio are against the adjusted-table run.

| Parameterisation | Step size | Leapfrogs | Divergences | Min ESS | Max Rhat | ESS per 1,000 gradients | Seconds | Gap | sd ratio |
|---|---|---|---|---|---|---|---|---|---|
| between bounds | 0.00936 | 498 | 5 | 369 | 1.002 | 0.12 | 43 | 0.05 | 0.95-1.04 |
| cell logit | 0.0337 | 124 | 0 | 200 | 1.008 | 0.27 | 10 | 0.07 | 0.97-1.17 |
| adjusted row | 0.137 | 41 | 1 | 533 | 1.003 | 2.17 | 8 | 0.08 | 0.90-1.03 |
| adjusted table | 0.142 | 34 | 0 | 786 | 1.001 | 3.82 | 30 | 0.00 | 1.00-1.00 |

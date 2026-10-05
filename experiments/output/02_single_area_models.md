# The single-area model in six coordinate systems

One table, placeholder prior (ILR coordinates ~ normal(0, 2)), short runs. Times are CPU seconds in the sampling phase, summed over the two chains, on a 2-core cloud machine; treat them as relative.

## Part 1: the stated Jacobian against finite differences

Six random points per row. 'Offset' is the mean of (finite difference - stated); it should be zero, or a constant for the ILR parameterisation, whose Jacobian is stated up to a constant. 'Spread' is the largest deviation from that mean and should be near zero.

| Size | Margin scaling | Parameterisation | Offset | Spread |
|---|---|---|---|---|
| 3x3 | off | ILR + log volume | 1.10e+00 | 9.3e-11 |
| 3x3 | off | log cells | -1.87e-10 | 1.7e-11 |
| 3x3 | off | between bounds | 1.46e-11 | 5.9e-11 |
| 3x3 | off | cell logit | 4.95e-11 | 9.7e-11 |
| 3x3 | off | adjusted row | 4.50e-11 | 4.4e-11 |
| 3x3 | off | adjusted table | 2.07e-11 | 1.6e-10 |
| 3x3 | on | between bounds | -2.34e-09 | 8.3e-09 |
| 3x3 | on | cell logit | -2.36e-09 | 6.0e-09 |
| 3x3 | on | adjusted row | 3.29e-09 | 6.9e-09 |
| 3x3 | on | adjusted table | 1.15e-09 | 4.1e-08 |
| 7x7 | off | ILR + log volume | 1.95e+00 | 2.6e-10 |
| 7x7 | off | log cells | -3.44e-10 | 9.6e-11 |
| 7x7 | off | between bounds | -1.39e-10 | 1.8e-10 |
| 7x7 | off | cell logit | 1.94e-10 | 3.1e-10 |
| 7x7 | off | adjusted row | 3.11e-10 | 2.1e-10 |
| 7x7 | off | adjusted table | 1.31e-10 | 8.7e-10 |
| 7x7 | on | between bounds | -4.89e-10 | 9.6e-09 |
| 7x7 | on | cell logit | 2.19e-09 | 5.7e-09 |
| 7x7 | on | adjusted row | -1.40e-10 | 8.5e-09 |
| 7x7 | on | adjusted table | 4.92e-08 | 6.7e-08 |

## Part 2: sampling the same posterior

2 chains, 500 warm-up + 500 draws each, started at the independence table, margin scaling off. Summaries are of the log cells; ESS is the smallest over cells, of 1,000 draws. 'Gap' is the largest difference in a posterior mean from the adjusted-row run, in posterior sds; 'sd ratio' is the range of the ratio of posterior sds to that run.

| Size | eps | Parameterisation | Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | ms per gradient | ESS per second | Gap | sd ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3x3 | 1 | ILR + log volume | 0.00845 | 420 | 0 | 1% | 144 | 1.007 | 0.34 | 0.003 | 112.0 | 0.12 | 0.92-1.17 |
| 3x3 | 1 | log cells | 0.00727 | 457 | 0 | 5% | 83 | 1.027 | 0.18 | 0.003 | 65.0 | 0.22 | 0.93-1.24 |
| 3x3 | 1 | between bounds | 0.191 | 22 | 0 | 0% | 264 | 1.012 | 12.02 | 0.008 | 1532.2 | 0.10 | 0.91-1.10 |
| 3x3 | 1 | cell logit | 0.276 | 13 | 0 | 0% | 9 | 1.251 | 0.72 | 0.009 | 82.5 | 0.46 | 0.93-1.61 |
| 3x3 | 1 | adjusted row | 0.351 | 11 | 0 | 0% | 503 | 1.002 | 47.39 | 0.021 | 2246.8 | 0.00 | 1.00-1.00 |
| 3x3 | 1 | adjusted table | 0.32 | 12 | 0 | 0% | 563 | 1.006 | 45.07 | 0.032 | 1427.9 | 0.16 | 0.95-1.05 |
| 3x3 | 0.1 | ILR + log volume | 0.000777 | 894 | 0 | 82% | 9 | 1.478 | 0.01 | 0.003 | 3.2 | 0.40 | 0.54-1.20 |
| 3x3 | 0.1 | log cells | 0.00163 | 790 | 83 | 67% | 3 | 1.254 | 0.00 | 0.002 | 1.5 | 1.14 | 0.23-1.18 |
| 3x3 | 0.1 | between bounds | 0.0941 | 34 | 0 | 0% | 335 | 1.008 | 9.75 | 0.007 | 1362.5 | 0.14 | 0.89-1.03 |
| 3x3 | 0.1 | cell logit | 0.131 | 21 | 0 | 0% | 154 | 1.006 | 7.39 | 0.008 | 919.8 | 0.10 | 0.91-1.17 |
| 3x3 | 0.1 | adjusted row | 0.145 | 19 | 0 | 0% | 215 | 1.016 | 11.62 | 0.014 | 854.9 | 0.00 | 1.00-1.00 |
| 3x3 | 0.1 | adjusted table | 0.151 | 18 | 0 | 0% | 260 | 1.022 | 14.73 | 0.026 | 560.8 | 0.13 | 0.93-1.19 |
| 5x5 | 1 | ILR + log volume | 0.00612 | 768 | 0 | 30% | 161 | 1.025 | 0.21 | 0.006 | 37.7 | 0.12 | 0.92-1.13 |
| 5x5 | 1 | log cells | 0.0108 | 468 | 10 | 0% | 68 | 1.031 | 0.14 | 0.005 | 29.9 | 0.36 | 0.38-1.09 |
| 5x5 | 1 | between bounds | 0.046 | 85 | 0 | 0% | 193 | 1.017 | 2.28 | 0.014 | 167.9 | 0.13 | 0.87-1.07 |
| 5x5 | 1 | cell logit | 0.136 | 30 | 0 | 0% | 181 | 1.019 | 5.96 | 0.013 | 445.0 | 0.10 | 0.91-1.08 |
| 5x5 | 1 | adjusted row | 0.219 | 19 | 0 | 0% | 216 | 1.003 | 11.25 | 0.028 | 397.6 | 0.00 | 1.00-1.00 |
| 5x5 | 1 | adjusted table | 0.215 | 25 | 0 | 0% | 298 | 1.003 | 11.96 | 0.066 | 180.2 | 0.15 | 0.88-1.12 |
| 5x5 | 0.1 | ILR + log volume | 0.00113 | 1023 | 0 | 100% | 2 | 1.682 | 0.00 | 0.006 | 0.3 | 0.27 | 0.75-1.27 |
| 5x5 | 0.1 | log cells | 0.00163 | 1022 | 1 | 100% | 1 | 1.663 | 0.00 | 0.004 | 0.3 | 1.01 | 0.17-1.42 |
| 5x5 | 0.1 | between bounds | 0.0228 | 173 | 0 | 0% | 180 | 1.015 | 1.04 | 0.011 | 96.7 | 0.18 | 0.76-1.04 |
| 5x5 | 0.1 | cell logit | 0.0933 | 42 | 0 | 0% | 121 | 1.018 | 2.88 | 0.017 | 166.6 | 0.20 | 0.79-1.12 |
| 5x5 | 0.1 | adjusted row | 0.129 | 30 | 0 | 0% | 165 | 1.007 | 5.50 | 0.025 | 219.2 | 0.00 | 1.00-1.00 |
| 5x5 | 0.1 | adjusted table | 0.137 | 28 | 0 | 0% | 235 | 1.021 | 8.46 | 0.069 | 123.3 | 0.16 | 0.94-1.07 |
| 7x7 | 1 | ILR + log volume | 0.0109 | 485 | 8 | 0% | 116 | 1.019 | 0.24 | 0.013 | 18.4 | 0.23 | 0.20-1.11 |
| 7x7 | 1 | log cells | 0.0112 | 474 | 0 | 0% | 61 | 1.036 | 0.13 | 0.008 | 16.0 | 0.28 | 0.13-1.12 |
| 7x7 | 1 | between bounds | 0.0164 | 255 | 0 | 0% | 112 | 1.017 | 0.44 | 0.021 | 21.1 | 0.22 | 0.87-1.09 |
| 7x7 | 1 | cell logit | 0.0752 | 57 | 0 | 0% | 19 | 1.118 | 0.33 | 0.019 | 17.6 | 0.35 | 0.90-1.56 |
| 7x7 | 1 | adjusted row | 0.133 | 43 | 0 | 0% | 124 | 1.004 | 2.89 | 0.042 | 69.1 | 0.00 | 1.00-1.00 |
| 7x7 | 1 | adjusted table | 0.178 | 30 | 0 | 0% | 127 | 1.021 | 4.17 | 0.196 | 21.2 | 0.28 | 0.83-1.08 |
| 7x7 | 0.1 | ILR + log volume | 0.00124 | 1023 | 0 | 100% | 20 | 1.131 | 0.02 | 0.012 | 1.6 | 0.54 | 0.14-1.29 |
| 7x7 | 0.1 | log cells | 0.00129 | 1023 | 0 | 100% | 17 | 1.082 | 0.02 | 0.008 | 2.1 | 0.40 | 0.04-1.14 |
| 7x7 | 0.1 | between bounds | 0.0107 | 414 | 0 | 0% | 33 | 1.048 | 0.08 | 0.020 | 4.0 | 0.20 | 0.90-1.20 |
| 7x7 | 0.1 | cell logit | 0.0629 | 64 | 0 | 0% | 6 | 1.195 | 0.09 | 0.022 | 4.3 | 0.26 | 0.93-1.34 |
| 7x7 | 0.1 | adjusted row | 0.12 | 43 | 0 | 0% | 66 | 1.016 | 1.53 | 0.040 | 38.4 | 0.00 | 1.00-1.00 |
| 7x7 | 0.1 | adjusted table | 0.125 | 48 | 0 | 0% | 155 | 1.024 | 3.25 | 0.181 | 18.0 | 0.12 | 0.89-1.15 |

## Part 3: longer runs of the sequential parameterisations at 7x7

2 chains, 1,000 warm-up + 3,000 draws each, eps = 1, margin scaling off. ESS is of 6,000 draws; gap and sd ratio are against the adjusted-table run.

| Parameterisation | Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | ms per gradient | ESS per second | Gap | sd ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| between bounds | 0.012 | 382 | 6 | 0% | 403 | 1.002 | 0.18 | 0.021 | 8.6 | 0.06 | 0.89-1.03 |
| cell logit | 0.0445 | 95 | 0 | 0% | 162 | 1.008 | 0.28 | 0.023 | 12.2 | 0.20 | 0.96-1.27 |
| adjusted row | 0.151 | 37 | 5 | 0% | 518 | 1.009 | 2.33 | 0.045 | 52.3 | 0.04 | 0.97-1.03 |
| adjusted table | 0.162 | 30 | 0 | 0% | 957 | 1.003 | 5.26 | 0.170 | 30.9 | 0.00 | 1.00-1.00 |

## Part 4: does scaling the margin parameters by the penalty width matter?

Adjusted row, 2 chains, 500 warm-up + 500 draws. 'off': the margin parameters are plain log deviations from the observed margins. 'on': they are divided by the penalty's width, so they are about unit scale.

| Size | eps | Margin scaling | Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | ms per gradient | ESS per second | Warm-up seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3x3 | 1 | off | 0.351 | 11 | 0 | 0% | 503 | 1.002 | 47.39 | 0.015 | 3145.6 | 1.4 |
| 3x3 | 1 | on | 0.373 | 11 | 0 | 0% | 387 | 1.006 | 34.23 | 0.013 | 2548.6 | 0.2 |
| 3x3 | 0.1 | off | 0.145 | 19 | 0 | 0% | 215 | 1.016 | 11.62 | 0.013 | 872.2 | 2.0 |
| 3x3 | 0.1 | on | 0.367 | 12 | 0 | 0% | 514 | 1.001 | 44.17 | 0.013 | 3520.2 | 0.1 |
| 3x3 | 0.01 | off | 0.0145 | 146 | 0 | 0% | 154 | 1.014 | 1.05 | 0.011 | 94.3 | 4.9 |
| 3x3 | 0.01 | on | 0.345 | 12 | 0 | 0% | 538 | 1.000 | 43.91 | 0.012 | 3662.2 | 0.1 |
| 7x7 | 1 | off | 0.133 | 43 | 0 | 0% | 124 | 1.004 | 2.89 | 0.039 | 74.4 | 7.8 |
| 7x7 | 1 | on | 0.17 | 31 | 0 | 0% | 54 | 1.030 | 1.76 | 0.043 | 41.0 | 1.6 |
| 7x7 | 0.1 | off | 0.12 | 43 | 0 | 0% | 66 | 1.016 | 1.53 | 0.042 | 36.8 | 9.2 |
| 7x7 | 0.1 | on | 0.171 | 30 | 0 | 0% | 65 | 1.037 | 2.15 | 0.042 | 51.7 | 1.4 |
| 7x7 | 0.01 | off | 0.0142 | 309 | 2 | 0% | 85 | 1.015 | 0.28 | 0.037 | 7.4 | 23.9 |
| 7x7 | 0.01 | on | 0.171 | 29 | 1 | 0% | 63 | 1.023 | 2.17 | 0.041 | 52.9 | 1.3 |

## Part 5: cost of one gradient by table size

Milliseconds per gradient, mean of the two Part 2 runs. The last column is the 7x7 cost over the 3x3 cost; the number of parameters grows 5.4 times.

| Parameterisation | 3x3 | 5x5 | 7x7 | 7x7 / 3x3 |
|---|---|---|---|---|
| ILR + log volume | 0.003 | 0.006 | 0.012 | 3.9 |
| log cells | 0.003 | 0.005 | 0.008 | 3.1 |
| between bounds | 0.008 | 0.012 | 0.020 | 2.7 |
| cell logit | 0.008 | 0.015 | 0.020 | 2.4 |
| adjusted row | 0.017 | 0.027 | 0.041 | 2.4 |
| adjusted table | 0.029 | 0.067 | 0.189 | 6.5 |

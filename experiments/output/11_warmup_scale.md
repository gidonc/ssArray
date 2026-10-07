# Warm-up and margin scale, Scottish grand total (eps 0.01, sigma_b 1.0)

|                                      |   failed |   max_rhat |   min_ess |   leapfrogs |   seconds |   metric_off_by |
|:-------------------------------------|---------:|-----------:|----------:|------------:|----------:|----------------:|
| ('position logit', 'default')        |        2 |       3.25 |      1.17 |     1023    |     24.73 |          402.5  |
| ('position logit', 'long warm-up')   |        0 |       1.01 |    808.54 |      988.18 |     62.69 |          173.37 |
| ('position logit', 'init metric')    |        0 |       1.03 |    129.25 |     1023    |     23.31 |          402.5  |
| ('position logit', 'scaled')         |        0 |       1.01 |    526.69 |       15.44 |      0.75 |            1.17 |
| ('odds-ratio logit', 'default')      |        2 |       2.12 |      1.55 |     1023    |     23.55 |          402.5  |
| ('odds-ratio logit', 'long warm-up') |        0 |       1.01 |    701.83 |     1010.2  |     58.04 |          173.37 |
| ('odds-ratio logit', 'init metric')  |        0 |       1.03 |    145.35 |     1023    |     21.33 |          402.5  |
| ('odds-ratio logit', 'scaled')       |        0 |       1    |    780.59 |       15.14 |      0.76 |            1.28 |

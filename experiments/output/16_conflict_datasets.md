# Conflict points on two elections (5 areas x 2 seeds per cell)

## Median ESS per 1000 gradients, converged runs

|                                                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('Scotland 2007', 0.5, 'indep', 'equal sd')                |                6.7 |         6   |             78.5 |               95.2 |           44.6 |             25.2 |
| ('Scotland 2007', 0.5, 'raked', 'equal sd')                |                3.5 |         8.5 |             67.4 |               78   |           42.2 |             17.9 |
| ('Scotland 2007', 0.5, 'national', 'equal sd')             |                2.6 |         7.1 |             63.2 |               57.8 |           37.2 |             23.8 |
| ('Scotland 2007', 0.5, 'national', 'interactions only')    |                1.8 |         4.7 |             60.7 |               57.8 |           42.4 |             24.8 |
| ('Scotland 2007', 2.0, 'indep', 'equal sd')                |                0.4 |         0.4 |              7.4 |               12.5 |           24.1 |             13.7 |
| ('Scotland 2007', 2.0, 'raked', 'equal sd')                |                0.4 |         0.1 |              3.2 |                8.1 |           12   |             11.8 |
| ('Scotland 2007', 2.0, 'national', 'equal sd')             |                0.3 |         0.2 |              2.7 |                6.4 |           10.5 |             12.9 |
| ('Scotland 2007', 2.0, 'national', 'interactions only')    |                0.3 |         0.1 |              1.1 |                2.9 |           13.8 |             19.4 |
| ('New Zealand 2017', 0.5, 'indep', 'equal sd')             |                3.9 |         4   |             53.9 |               41.7 |           30.2 |             14   |
| ('New Zealand 2017', 0.5, 'raked', 'equal sd')             |                3.2 |         8.5 |             37.9 |               41.3 |           24.9 |             13.6 |
| ('New Zealand 2017', 0.5, 'national', 'equal sd')          |                3.2 |         7.4 |             45.3 |               38.6 |           28.2 |             16   |
| ('New Zealand 2017', 0.5, 'national', 'interactions only') |                1.4 |         4.5 |             29   |               29.6 |           26.8 |             20   |
| ('New Zealand 2017', 2.0, 'indep', 'equal sd')             |                0.4 |         0.2 |              4.3 |                7.8 |           21.6 |             11.5 |
| ('New Zealand 2017', 2.0, 'raked', 'equal sd')             |                0.4 |         0.5 |              1.8 |                5.2 |           13.5 |             10.2 |
| ('New Zealand 2017', 2.0, 'national', 'equal sd')          |                0.5 |         0.3 |              2.5 |                5.5 |            9.2 |             11.6 |
| ('New Zealand 2017', 2.0, 'national', 'interactions only') |                0.2 |         0.3 |              0.6 |                2.4 |           11.3 |             15   |

## Runs with rhat > 1.05, of 10

|                                                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('Scotland 2007', 0.5, 'indep', 'equal sd')                |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 0.5, 'raked', 'equal sd')                |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 0.5, 'national', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 0.5, 'national', 'interactions only')    |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 2.0, 'indep', 'equal sd')                |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 2.0, 'raked', 'equal sd')                |                  0 |           3 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 2.0, 'national', 'equal sd')             |                  0 |           3 |                0 |                  0 |              0 |                0 |
| ('Scotland 2007', 2.0, 'national', 'interactions only')    |                  0 |           2 |                0 |                  1 |              1 |                0 |
| ('New Zealand 2017', 0.5, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'raked', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'national', 'equal sd')          |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'national', 'interactions only') |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'raked', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'national', 'equal sd')          |                  0 |           1 |                0 |                  0 |              1 |                0 |
| ('New Zealand 2017', 2.0, 'national', 'interactions only') |                  0 |           2 |                0 |                  0 |              0 |                0 |

## Volume (middle area, raked centre, sigma_b 1): median ESS per 1000 gradients

|                                 |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:--------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('New Zealand 2017', 354.22)    |              14.36 |       28.55 |             3.8  |              16.49 |          12.6  |             7.69 |
| ('New Zealand 2017', 35422.0)   |               1.25 |        2.18 |            22.58 |              28.98 |          21.62 |            12.95 |
| ('New Zealand 2017', 3542200.0) |             nan    |        0.25 |            21.84 |              26.35 |          21.45 |            11.2  |
| ('Scotland 2007', 296.92)       |              18.41 |       15.69 |            24.42 |              26.95 |          23.86 |            16.98 |
| ('Scotland 2007', 29692.0)      |               1.13 |        1.13 |            42.61 |              49.52 |          31.5  |            16.98 |
| ('Scotland 2007', 2969200.0)    |             nan    |        0.11 |            45.7  |              35.35 |          34.56 |            15.1  |

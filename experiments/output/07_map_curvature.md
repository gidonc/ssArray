# Nonlinearity of the coordinate maps over one posterior-sd step (random whitened directions)

r = |actual change - linear prediction| / |linear prediction|; 0 = linear. Rows with rhat >= 1.05 have an unreliable posterior covariance.
ilr = coordinates the prior is Gaussian in; margins = log row and column sums the penalty acts on.

|                                         |   ('median', 'ilr') |   ('median', 'margins') |   ('p90', 'ilr') |   ('p90', 'margins') |
|:----------------------------------------|--------------------:|------------------------:|-----------------:|---------------------:|
| (0.1, 1.0, 'ILR + log volume', 4.25606) |               0     |                   0.086 |            0     |                0.148 |
| (0.1, 1.0, 'adjusted row', 1.00325)     |               0.037 |                   0     |            0.055 |                0     |
| (0.1, 1.0, 'adjusted table', 1.00392)   |               0.041 |                   0     |            0.061 |                0     |
| (0.1, 1.0, 'log cells', 3.23011)        |               0     |                   0.118 |            0     |                0.201 |
| (0.1, 1.0, 'odds-ratio logit', 1.00146) |               0.049 |                   0     |            0.072 |                0     |
| (0.1, 1.0, 'position logit', 1.00861)   |               0.056 |                   0     |            0.123 |                0     |
| (0.1, 2.0, 'ILR + log volume', 7.42254) |               0     |                   0.09  |            0     |                0.164 |
| (0.1, 2.0, 'adjusted row', 1.00666)     |               0.064 |                   0     |            0.104 |                0     |
| (0.1, 2.0, 'adjusted table', 1.01533)   |               0.073 |                   0     |            0.115 |                0     |
| (0.1, 2.0, 'log cells', 3.80149)        |               0     |                   0.132 |            0     |                0.233 |
| (0.1, 2.0, 'odds-ratio logit', 1.00752) |               0.096 |                   0     |            0.244 |                0     |
| (0.1, 2.0, 'position logit', 1.0069)    |               0.12  |                   0     |            0.53  |                0     |
| (1.0, 1.0, 'ILR + log volume', 1.02454) |               0     |                   0.094 |            0     |                0.159 |
| (1.0, 1.0, 'adjusted row', 1.00483)     |               0.037 |                   0     |            0.056 |                0     |
| (1.0, 1.0, 'adjusted table', 1.0051)    |               0.041 |                   0     |            0.062 |                0     |
| (1.0, 1.0, 'log cells', 1.02001)        |               0     |                   0.151 |            0     |                0.231 |
| (1.0, 1.0, 'odds-ratio logit', 1.0028)  |               0.05  |                   0     |            0.071 |                0     |
| (1.0, 1.0, 'position logit', 1.0038)    |               0.057 |                   0     |            0.121 |                0     |
| (1.0, 2.0, 'ILR + log volume', 1.64923) |               0     |                   0.115 |            0     |                0.205 |
| (1.0, 2.0, 'adjusted row', 1.021)       |               0.065 |                   0     |            0.106 |                0     |
| (1.0, 2.0, 'adjusted table', 1.00941)   |               0.073 |                   0     |            0.118 |                0     |
| (1.0, 2.0, 'log cells', 1.37402)        |               0     |                   0.166 |            0     |                0.273 |
| (1.0, 2.0, 'odds-ratio logit', 1.01368) |               0.092 |                   0     |            0.223 |                0     |
| (1.0, 2.0, 'position logit', 1.00761)   |               0.125 |                   0     |            0.572 |                0     |

Step 2 sd: see the csv.

# Size and margin dials: median ESS per 1000 gradients (two replicates)

| design                       |   total |   row_margin_mean |   cell_mean |   margin_min |   margin_max |   ILR + log volume |   adjusted row |   log cells |   odds-ratio logit |
|:-----------------------------|--------:|------------------:|------------:|-------------:|-------------:|-------------------:|---------------:|------------:|-------------------:|
| 5x5 baseline                 |   30000 |              6000 |        1200 |         4962 |         6862 |                9.2 |           43.3 |         8.1 |              151.3 |
| 3x3, total fixed             |   30000 |             10000 |        3333 |         8151 |        11230 |                4.4 |          137   |         2.8 |              179.6 |
| 3x3, per-margin fixed        |   18000 |              6000 |        2000 |         4891 |         6738 |                6.5 |          174.9 |         4.2 |              190.3 |
| 3x3, per-cell fixed          |   10800 |              3600 |        1200 |         2934 |         4043 |                7.8 |          166.3 |         5.3 |              187.9 |
| 8x8, total fixed             |   30000 |              3750 |         469 |         3064 |         4231 |               13.4 |           21.8 |         9.3 |               84.5 |
| 8x8, per-margin fixed        |   48000 |              6000 |         750 |         4903 |         6769 |               14.6 |           18.9 |        11.3 |               78   |
| 8x8, per-cell fixed          |   76800 |              9600 |        1200 |         7845 |        10830 |                9.1 |           21.1 |         8.4 |               63   |
| rows geometric 0.5           |   30000 |              6000 |        1200 |          737 |        17245 |                5.5 |           48.5 |         5.2 |              113.7 |
| rows one dominant 0.5        |   30000 |              6000 |        1200 |         2453 |        18211 |                5.5 |           43.4 |         4.1 |              122.4 |
| rows one tiny 0.85           |   30000 |              6000 |        1200 |         1132 |         7850 |                8.8 |           40.9 |         8.5 |              128.6 |
| both geometric 0.5, aligned  |   30000 |              6000 |        1200 |          616 |        17625 |                4.5 |           24   |         4.2 |               25.5 |
| both geometric 0.5, reversed |   30000 |              6000 |        1200 |          700 |        17282 |                6.3 |           46.2 |         3.7 |              143   |

## Leapfrogs per draw

| design                       |   ILR + log volume |   adjusted row |   log cells |   odds-ratio logit |
|:-----------------------------|-------------------:|---------------:|------------:|-------------------:|
| 5x5 baseline                 |                 54 |             10 |         125 |                  7 |
| 3x3, total fixed             |                 48 |              7 |         123 |                  7 |
| 3x3, per-margin fixed        |                 38 |              7 |          91 |                  7 |
| 3x3, per-cell fixed          |                 31 |              7 |          77 |                  7 |
| 8x8, total fixed             |                 38 |             12 |          82 |                 11 |
| 8x8, per-margin fixed        |                 61 |             15 |         123 |                 12 |
| 8x8, per-cell fixed          |                 63 |             13 |         127 |                 12 |
| rows geometric 0.5           |                124 |              9 |         121 |                  8 |
| rows one dominant 0.5        |                122 |              9 |          85 |                  7 |
| rows one tiny 0.85           |                 67 |             10 |         125 |                  8 |
| both geometric 0.5, aligned  |                148 |             21 |         101 |                 19 |
| both geometric 0.5, reversed |                126 |             10 |         105 |                  7 |

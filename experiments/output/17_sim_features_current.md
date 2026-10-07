# Simulated designs, one feature at a time, current model: median over three replicates

Direct schemes from the earlier run.

## ESS per 1000 gradients

| arm                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| baseline                       |                6.4 |         5.4 |             92   |              112.2 |          171.8 |            137.4 |
| low volume                     |              172   |       117.7 |             47.7 |               61.2 |           58.3 |             67.2 |
| high volume                    |                0.6 |         0.6 |             87   |               91.4 |          149.8 |            168.4 |
| wide prior                     |                0.5 |         0.5 |              5.2 |                5.5 |           49.1 |             56.8 |
| diagonal centre                |                3   |         9   |             67.2 |               57.7 |          129.5 |            141.2 |
| diagonal centre, wide prior    |                0.3 |         0.3 |              2.2 |                5   |           22.7 |             41.1 |
| tied margins, diagonal, wide   |                0.3 |         2.1 |              1   |                3.3 |           15   |             33   |
| tied margins, wide             |                0.4 |         0.4 |              3   |                9.6 |           51.6 |             67   |
| tiny margin                    |                6.5 |         4.5 |             56.6 |               87.2 |          125.9 |            169.7 |
| tiny margin, interactions-only |                1.1 |         0.7 |              3.5 |                2.4 |            2.9 |            nan   |
| interactions-only, wide        |                0.3 |         0.2 |              2.3 |                6.3 |           45.7 |            117.8 |
| conflict 5 sd                  |                6.4 |         4.3 |             62.8 |               53.3 |          128.1 |            150.7 |
| conflict 20 sd                 |                0.4 |         2.5 |             26.9 |               20.5 |           19.6 |             45.5 |
| 3x3                            |                3.6 |         2.6 |            181   |              137.5 |          170.4 |            181.5 |
| 8x8                            |                9.6 |         4.7 |             66.9 |               54.2 |           78.9 |             91.2 |

## ESS per second

| arm                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| baseline                       |                307 |         527 |             1863 |               2762 |           2766 |             1368 |
| low volume                     |               4887 |        3847 |             1682 |               2150 |           1187 |             1002 |
| high volume                    |                 66 |          83 |             2139 |               1858 |           2684 |             2001 |
| wide prior                     |                 50 |          60 |              215 |                218 |            922 |              617 |
| diagonal centre                |                252 |         523 |             1775 |               1588 |           1647 |             1018 |
| diagonal centre, wide prior    |                 37 |          27 |              101 |                181 |            434 |              351 |
| tied margins, diagonal, wide   |                 33 |          70 |               41 |                126 |            209 |              274 |
| tied margins, wide             |                 42 |          58 |              154 |                380 |           1103 |              722 |
| tiny margin                    |                472 |         484 |             1966 |               2092 |           2209 |             1763 |
| tiny margin, interactions-only |                101 |          70 |               90 |                108 |             43 |              nan |
| interactions-only, wide        |                 28 |          33 |              107 |                174 |            833 |             1077 |
| conflict 5 sd                  |                421 |         489 |             1758 |               1633 |           2121 |             1579 |
| conflict 20 sd                 |                 29 |         193 |              734 |                623 |            344 |              453 |
| 3x3                            |                405 |         454 |             6554 |               4525 |           6304 |             3594 |
| 8x8                            |                293 |         196 |             1143 |               1041 |            812 |              534 |

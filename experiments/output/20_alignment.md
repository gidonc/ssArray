# Alignment and the left-over cell: median ESS per 1000 gradients (three replicates)

|                                                                |   left_over_share |   last_row_share |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:---------------------------------------------------------------|------------------:|-----------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('aligned, balance 1.0', 'smallest first', 'diagonal')         |             0.185 |            0.244 |             72.5 |              142.4 |           44.5 |             24   |
| ('aligned, balance 1.0', 'smallest first', 'largest column')   |             0.218 |            0.244 |             84.3 |              120.3 |           46.5 |            nan   |
| ('aligned, balance 1.0', 'smallest first', 'smallest column')  |             0.184 |            0.244 |             63.3 |              139.3 |           46.2 |            nan   |
| ('aligned, balance 0.75', 'smallest first', 'diagonal')        |             0.072 |            0.411 |             71.5 |               69.5 |           47.7 |             18.6 |
| ('aligned, balance 0.75', 'smallest first', 'largest column')  |             0.404 |            0.411 |             68.3 |               70.6 |           52.4 |            nan   |
| ('aligned, balance 0.75', 'smallest first', 'smallest column') |             0.072 |            0.411 |             65.5 |               79.2 |           35.7 |            nan   |
| ('aligned, balance 0.5', 'smallest first', 'diagonal')         |             0.022 |            0.586 |             46.3 |               25.8 |           25   |             16.5 |
| ('aligned, balance 0.5', 'smallest first', 'largest column')   |             0.588 |            0.586 |             43.8 |               49.1 |           36.9 |            nan   |
| ('aligned, balance 0.5', 'smallest first', 'smallest column')  |             0.022 |            0.586 |             46.8 |               43.1 |           23   |            nan   |
| ('aligned, balance 0.5', 'largest first', 'diagonal')          |             0.048 |            0.023 |             20.2 |               30.6 |           24.7 |            nan   |
| ('aligned, balance 0.5', 'largest first', 'largest column')    |             0.588 |            0.023 |             45.8 |               21.3 |           25.5 |            nan   |
| ('aligned, balance 0.5', 'largest first', 'smallest column')   |             0.022 |            0.023 |             19.1 |               34.6 |           20.4 |            nan   |
| ('aligned, balance 0.35', 'smallest first', 'diagonal')        |             0.004 |            0.742 |             17.9 |               13.9 |           13.3 |              7.6 |
| ('aligned, balance 0.35', 'smallest first', 'largest column')  |             0.75  |            0.742 |             19.5 |               15.8 |           14.7 |            nan   |
| ('aligned, balance 0.35', 'smallest first', 'smallest column') |             0.004 |            0.742 |             19.4 |               15.3 |           11.6 |            nan   |
| ('reversed, balance 0.5', 'smallest first', 'diagonal')        |             0.045 |            0.586 |            148.2 |              136.4 |           47.4 |             29.1 |
| ('reversed, balance 0.5', 'smallest first', 'largest column')  |             0.567 |            0.586 |            162.2 |              136   |           69.6 |            nan   |
| ('reversed, balance 0.5', 'smallest first', 'smallest column') |             0.022 |            0.586 |            131   |              130.5 |           38.4 |            nan   |
| ('shifted, balance 0.5', 'smallest first', 'diagonal')         |             0.051 |            0.586 |             72.8 |               63   |           39.1 |             15.6 |
| ('shifted, balance 0.5', 'smallest first', 'largest column')   |             0.542 |            0.586 |             74.8 |               58.6 |           43.9 |            nan   |
| ('shifted, balance 0.5', 'smallest first', 'smallest column')  |             0.024 |            0.586 |             63.4 |               73.9 |           36.6 |            nan   |

## Leapfrogs per draw

|                                                                |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:---------------------------------------------------------------|-----------------:|-------------------:|---------------:|-----------------:|
| ('aligned, balance 1.0', 'smallest first', 'diagonal')         |               10 |                  7 |             10 |               14 |
| ('aligned, balance 1.0', 'smallest first', 'largest column')   |                8 |                  7 |              9 |              nan |
| ('aligned, balance 1.0', 'smallest first', 'smallest column')  |               11 |                  7 |             10 |              nan |
| ('aligned, balance 0.75', 'smallest first', 'diagonal')        |               14 |                 11 |             13 |               15 |
| ('aligned, balance 0.75', 'smallest first', 'largest column')  |               13 |                 13 |             11 |              nan |
| ('aligned, balance 0.75', 'smallest first', 'smallest column') |               15 |                 12 |             12 |              nan |
| ('aligned, balance 0.5', 'smallest first', 'diagonal')         |               15 |                 18 |             21 |               25 |
| ('aligned, balance 0.5', 'smallest first', 'largest column')   |               15 |                 15 |             16 |              nan |
| ('aligned, balance 0.5', 'smallest first', 'smallest column')  |               16 |                 17 |             19 |              nan |
| ('aligned, balance 0.5', 'largest first', 'diagonal')          |               25 |                 22 |             23 |              nan |
| ('aligned, balance 0.5', 'largest first', 'largest column')    |               17 |                 19 |             17 |              nan |
| ('aligned, balance 0.5', 'largest first', 'smallest column')   |               27 |                 19 |             19 |              nan |
| ('aligned, balance 0.35', 'smallest first', 'diagonal')        |               35 |                 34 |             37 |               38 |
| ('aligned, balance 0.35', 'smallest first', 'largest column')  |               32 |                 34 |             35 |              nan |
| ('aligned, balance 0.35', 'smallest first', 'smallest column') |               34 |                 48 |             40 |              nan |
| ('reversed, balance 0.5', 'smallest first', 'diagonal')        |                7 |                  7 |             10 |               14 |
| ('reversed, balance 0.5', 'smallest first', 'largest column')  |                7 |                  7 |              8 |              nan |
| ('reversed, balance 0.5', 'smallest first', 'smallest column') |                7 |                  7 |             11 |              nan |
| ('shifted, balance 0.5', 'smallest first', 'diagonal')         |               14 |                 13 |             16 |               16 |
| ('shifted, balance 0.5', 'smallest first', 'largest column')   |               14 |                 14 |             14 |              nan |
| ('shifted, balance 0.5', 'smallest first', 'smallest column')  |               14 |                 15 |             16 |              nan |

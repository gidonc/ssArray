# Four datasets with the current model (largest-column margin coordinates, centred interior, split basis); five areas x two seeds

Direct schemes (ILR, log cells) are taken from the earlier run: the model changes do not touch them.

## Median ESS per 1000 gradients (converged runs)

|                                                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('Scotland 2007', 0.5, 'indep', 'equal sd')                |                6.7 |         6   |            119.5 |              124.3 |          171.4 |            170   |
| ('Scotland 2007', 0.5, 'raked', 'equal sd')                |                3.5 |         8.5 |             84.1 |              102.8 |          136.5 |            121.9 |
| ('Scotland 2007', 0.5, 'national', 'equal sd')             |                2.6 |         7.1 |             83.3 |               74.5 |           98.9 |             78.2 |
| ('Scotland 2007', 0.5, 'national', 'interactions only')    |                1.8 |         4.7 |             69.3 |               72.4 |           90.9 |            188.1 |
| ('Scotland 2007', 2.0, 'indep', 'equal sd')                |                0.4 |         0.4 |              6.3 |               13.2 |           50.3 |             50.8 |
| ('Scotland 2007', 2.0, 'raked', 'equal sd')                |                0.4 |         0.1 |              2.9 |                6.7 |           24.9 |             37.6 |
| ('Scotland 2007', 2.0, 'national', 'equal sd')             |                0.3 |         0.2 |              2.6 |                5.6 |           24.2 |             35   |
| ('Scotland 2007', 2.0, 'national', 'interactions only')    |                0.3 |         0.1 |              1.2 |                3   |           32.2 |             81.6 |
| ('New Zealand 2017', 0.5, 'indep', 'equal sd')             |                3.9 |         4   |             89.4 |               55.8 |          135.6 |            139   |
| ('New Zealand 2017', 0.5, 'raked', 'equal sd')             |                3.2 |         8.5 |             61.3 |               57.3 |          121   |            135.7 |
| ('New Zealand 2017', 0.5, 'national', 'equal sd')          |                3.2 |         7.4 |             59.6 |               58.6 |          123.3 |            125.3 |
| ('New Zealand 2017', 0.5, 'national', 'interactions only') |                1.4 |         4.5 |             49.6 |               43.6 |          104.9 |            190   |
| ('New Zealand 2017', 2.0, 'indep', 'equal sd')             |                0.4 |         0.2 |              6.4 |                9.6 |           47.4 |             51.6 |
| ('New Zealand 2017', 2.0, 'raked', 'equal sd')             |                0.4 |         0.5 |              2.2 |                4.7 |           36.1 |             41.4 |
| ('New Zealand 2017', 2.0, 'national', 'equal sd')          |                0.5 |         0.3 |              2.4 |                6.6 |           33.9 |             44   |
| ('New Zealand 2017', 2.0, 'national', 'interactions only') |                0.2 |         0.3 |              0.6 |                3.3 |           32   |             93.7 |
| ('senc', 0.5, 'indep', 'equal sd')                         |               13.6 |        15.6 |            171.2 |              155.9 |          162   |            166.5 |
| ('senc', 0.5, 'raked', 'equal sd')                         |               14.1 |        46   |            156.4 |              157.4 |          152.5 |            134.4 |
| ('senc', 0.5, 'national', 'equal sd')                      |               14.2 |        65.1 |            156.8 |              157.1 |          190.3 |            158.8 |
| ('senc', 0.5, 'national', 'interactions only')             |                2.8 |         0.7 |             11.2 |               12.1 |            7   |            175.1 |
| ('senc', 2.0, 'indep', 'equal sd')                         |                1.8 |         1.7 |             72.5 |               52.7 |           82.8 |             76.6 |
| ('senc', 2.0, 'raked', 'equal sd')                         |                2   |         1.7 |             45.9 |               47   |           56   |             63   |
| ('senc', 2.0, 'national', 'equal sd')                      |                1.8 |         1.6 |             55.6 |               45.5 |           55.7 |             71.4 |
| ('senc', 2.0, 'national', 'interactions only')             |                0.7 |         1.4 |              4.8 |                3.8 |            6.6 |             81.2 |
| ('redistrict', 0.5, 'indep', 'equal sd')                   |               26   |        22.1 |            175   |              178.2 |          184.6 |            177   |
| ('redistrict', 0.5, 'national', 'equal sd')                |               26.2 |        20.1 |            180   |              164.7 |          143.3 |            146.3 |
| ('redistrict', 0.5, 'national', 'interactions only')       |               10.6 |         1.1 |            132.9 |               86.7 |           13.7 |             93.7 |
| ('redistrict', 2.0, 'indep', 'equal sd')                   |                3   |         1.9 |             40.7 |               32.2 |           71.9 |             78.5 |
| ('redistrict', 2.0, 'national', 'equal sd')                |                2.2 |         1.9 |             16.1 |               21.4 |           34.4 |             63.8 |
| ('redistrict', 2.0, 'national', 'interactions only')       |                0.5 |         0.9 |              5   |                8.9 |            8.3 |             95.6 |

## Median ESS per second

|                                                            |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('Scotland 2007', 0.5, 'indep', 'equal sd')                |                742 |        1037 |             2814 |               2890 |           2883 |             1698 |
| ('Scotland 2007', 0.5, 'raked', 'equal sd')                |                506 |        1016 |             2321 |               2614 |           1831 |             1073 |
| ('Scotland 2007', 0.5, 'national', 'equal sd')             |                391 |         839 |             1963 |               1862 |           1695 |              739 |
| ('Scotland 2007', 0.5, 'national', 'interactions only')    |                248 |         609 |             1853 |               2095 |           1539 |             1635 |
| ('Scotland 2007', 2.0, 'indep', 'equal sd')                |                 72 |          89 |              315 |                502 |            927 |              608 |
| ('Scotland 2007', 2.0, 'raked', 'equal sd')                |                 73 |          21 |              141 |                299 |            488 |              406 |
| ('Scotland 2007', 2.0, 'national', 'equal sd')             |                 60 |          41 |              119 |                230 |            441 |              398 |
| ('Scotland 2007', 2.0, 'national', 'interactions only')    |                 55 |          28 |               48 |                118 |            602 |              649 |
| ('New Zealand 2017', 0.5, 'indep', 'equal sd')             |                503 |         646 |             2175 |               1839 |           2093 |             1458 |
| ('New Zealand 2017', 0.5, 'raked', 'equal sd')             |                483 |         949 |             1891 |               1746 |           1828 |             1237 |
| ('New Zealand 2017', 0.5, 'national', 'equal sd')          |                460 |         815 |             2057 |               1843 |           1966 |             1093 |
| ('New Zealand 2017', 0.5, 'national', 'interactions only') |                198 |         502 |             1580 |               1455 |           1686 |             1873 |
| ('New Zealand 2017', 2.0, 'indep', 'equal sd')             |                 71 |          51 |              270 |                344 |           1028 |              640 |
| ('New Zealand 2017', 2.0, 'raked', 'equal sd')             |                 80 |          93 |               94 |                206 |            772 |              496 |
| ('New Zealand 2017', 2.0, 'national', 'equal sd')          |                 86 |          58 |              107 |                292 |            745 |              432 |
| ('New Zealand 2017', 2.0, 'national', 'interactions only') |                 45 |          61 |               31 |                139 |            641 |              731 |
| ('senc', 0.5, 'indep', 'equal sd')                         |               2521 |        3206 |             7140 |               6415 |           5102 |             3207 |
| ('senc', 0.5, 'raked', 'equal sd')                         |               2620 |        4378 |             6170 |               6192 |           4487 |             2453 |
| ('senc', 0.5, 'national', 'equal sd')                      |               2432 |        7031 |             6232 |               6200 |           5733 |             2634 |
| ('senc', 0.5, 'national', 'interactions only')             |                286 |         247 |              328 |                271 |            212 |             3367 |
| ('senc', 2.0, 'indep', 'equal sd')                         |                461 |         515 |             2946 |               2680 |           2733 |             1617 |
| ('senc', 2.0, 'raked', 'equal sd')                         |                568 |         324 |             2268 |               2296 |           1970 |             1127 |
| ('senc', 2.0, 'national', 'equal sd')                      |                485 |         322 |             1901 |               1659 |           1910 |             1242 |
| ('senc', 2.0, 'national', 'interactions only')             |                176 |         173 |              205 |                200 |            136 |             1349 |
| ('redistrict', 0.5, 'indep', 'equal sd')                   |               3930 |        4146 |             7951 |               7309 |           5787 |             3520 |
| ('redistrict', 0.5, 'national', 'equal sd')                |               3121 |        3703 |             6942 |               7351 |           4354 |             2847 |
| ('redistrict', 0.5, 'national', 'interactions only')       |               1556 |         401 |             5721 |               3266 |            402 |             1556 |
| ('redistrict', 2.0, 'indep', 'equal sd')                   |                717 |         637 |             1883 |               1569 |           2080 |             1452 |
| ('redistrict', 2.0, 'national', 'equal sd')                |                514 |         555 |              954 |               1348 |           1220 |              996 |
| ('redistrict', 2.0, 'national', 'interactions only')       |                119 |         311 |              234 |                369 |            268 |             1306 |

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
| ('Scotland 2007', 2.0, 'national', 'interactions only')    |                  0 |           2 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'raked', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'national', 'equal sd')          |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 0.5, 'national', 'interactions only') |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'raked', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'national', 'equal sd')          |                  0 |           1 |                0 |                  1 |              0 |                0 |
| ('New Zealand 2017', 2.0, 'national', 'interactions only') |                  0 |           2 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'indep', 'equal sd')                         |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'raked', 'equal sd')                         |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'equal sd')                      |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'interactions only')             |                  3 |           1 |                4 |                  5 |              2 |                8 |
| ('senc', 2.0, 'indep', 'equal sd')                         |                  0 |           1 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'raked', 'equal sd')                         |                  0 |           4 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'national', 'equal sd')                      |                  0 |           6 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'national', 'interactions only')             |                  1 |           9 |                3 |                  3 |              5 |                7 |
| ('redistrict', 0.5, 'indep', 'equal sd')                   |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('redistrict', 0.5, 'national', 'equal sd')                |                  2 |           2 |                2 |                  2 |              0 |                1 |
| ('redistrict', 0.5, 'national', 'interactions only')       |                  2 |           0 |                3 |                  2 |              1 |                2 |
| ('redistrict', 2.0, 'indep', 'equal sd')                   |                  0 |           1 |                0 |                  1 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'equal sd')                |                  0 |           3 |                0 |                  0 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'interactions only')       |                  0 |           4 |                2 |                  1 |              2 |                3 |

## Volume (middle area, sigma_b 1): median ESS per 1000 gradients

|                                   |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:----------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('New Zealand 2017', 354.22)      |               14.4 |        28.6 |             16.7 |               22.7 |           64.4 |             73.1 |
| ('New Zealand 2017', 35422.0)     |                1.3 |         2.2 |             35.1 |               21.4 |           89.4 |             77.2 |
| ('New Zealand 2017', 3542200.0)   |              nan   |         0.3 |             23.9 |               38.5 |           81.3 |             90.7 |
| ('Scotland 2007', 296.92)         |               18.4 |        15.7 |             21.7 |               29.6 |           49.6 |             59.1 |
| ('Scotland 2007', 29692.0)        |                1.1 |         1.1 |             40.7 |               37.3 |           64.2 |             95.6 |
| ('Scotland 2007', 2969200.0)      |              nan   |         0.1 |             46.1 |               48   |           74.5 |             73.4 |
| ('redistrict', 5.23)              |              103.8 |        63.5 |             60.8 |               35   |           44.6 |             36.4 |
| ('redistrict', 523.0000000000001) |                7.2 |        39.1 |            200.6 |              142.6 |          189.1 |            161.2 |
| ('redistrict', 52300.0)           |                0.8 |         2.5 |            174.4 |              152.1 |          165.9 |            163.8 |
| ('senc', 15.9)                    |               54.7 |        43.9 |             54.3 |               46.3 |           38.1 |             39.1 |
| ('senc', 1590.0)                  |                6.1 |         3.9 |             96   |               74.4 |           81.4 |             94.6 |
| ('senc', 159000.0)                |                0.6 |         0.7 |             89.3 |               80.4 |          110.1 |            106.2 |

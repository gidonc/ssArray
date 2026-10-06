# Conflict points on the eiPack data: senc (3x3, known cells) and redistrict (3x3, margins only; common centre = pooled independence)

## Median ESS per 1000 gradients, converged runs

|                                                      |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('senc', 0.5, 'indep', 'equal sd')                   |               13.6 |        15.6 |            130.1 |              125.8 |           99.4 |             61.9 |
| ('senc', 0.5, 'raked', 'equal sd')                   |               14.1 |        46   |            105.1 |               93   |           73.6 |             54.8 |
| ('senc', 0.5, 'national', 'equal sd')                |               14.2 |        65.1 |             97.2 |               96.2 |           70   |             51.2 |
| ('senc', 0.5, 'national', 'interactions only')       |                2.8 |         0.7 |              5.9 |                3.1 |            9.3 |             59   |
| ('senc', 2.0, 'indep', 'equal sd')                   |                1.8 |         1.7 |             45.1 |               35   |           56.4 |             44.7 |
| ('senc', 2.0, 'raked', 'equal sd')                   |                2   |         1.7 |             38   |               30.4 |           44.5 |             38.3 |
| ('senc', 2.0, 'national', 'equal sd')                |                1.8 |         1.6 |             52.4 |               41.1 |           55.6 |             36.5 |
| ('senc', 2.0, 'national', 'interactions only')       |                0.7 |         1.4 |              1.9 |               16.3 |            5.6 |             54.3 |
| ('redistrict', 0.5, 'indep', 'equal sd')             |               26   |        22.1 |            167.3 |              166.2 |          130.4 |             57.4 |
| ('redistrict', 0.5, 'national', 'equal sd')          |               26.2 |        20.1 |            164.8 |              181.6 |          104.2 |             48.7 |
| ('redistrict', 0.5, 'national', 'interactions only') |               10.6 |         1.1 |             19.9 |               87.5 |          152.4 |             72.2 |
| ('redistrict', 2.0, 'indep', 'equal sd')             |                3   |         1.9 |             40.9 |               35.9 |           33.6 |             37.7 |
| ('redistrict', 2.0, 'national', 'equal sd')          |                2.2 |         1.9 |             14.2 |               23.4 |           35   |             36.3 |
| ('redistrict', 2.0, 'national', 'interactions only') |                0.5 |         0.9 |              5.4 |                9.9 |           16.2 |             41.8 |

## Runs with rhat > 1.05 (of 10)

|                                                      |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('senc', 0.5, 'indep', 'equal sd')                   |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'raked', 'equal sd')                   |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'equal sd')                |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'interactions only')       |                  3 |           1 |                5 |                  2 |              4 |                7 |
| ('senc', 2.0, 'indep', 'equal sd')                   |                  0 |           1 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'raked', 'equal sd')                   |                  0 |           4 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'national', 'equal sd')                |                  0 |           6 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'national', 'interactions only')       |                  1 |           9 |                2 |                  8 |              4 |                7 |
| ('redistrict', 0.5, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('redistrict', 0.5, 'national', 'equal sd')          |                  2 |           2 |                2 |                  2 |              1 |                2 |
| ('redistrict', 0.5, 'national', 'interactions only') |                  2 |           0 |                1 |                  2 |              3 |                6 |
| ('redistrict', 2.0, 'indep', 'equal sd')             |                  0 |           1 |                0 |                  0 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'equal sd')          |                  0 |           3 |                0 |                  0 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'interactions only') |                  0 |           4 |                3 |                  2 |              1 |                5 |

## Divergences

|                                                      |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:-----------------------------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('senc', 0.5, 'indep', 'equal sd')                   |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'raked', 'equal sd')                   |                  0 |           1 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'equal sd')                |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('senc', 0.5, 'national', 'interactions only')       |               2247 |          51 |             3660 |               3113 |           3437 |             6117 |
| ('senc', 2.0, 'indep', 'equal sd')                   |                  0 |          22 |                0 |                  0 |              0 |                0 |
| ('senc', 2.0, 'raked', 'equal sd')                   |                  0 |         163 |                7 |                  0 |              0 |                0 |
| ('senc', 2.0, 'national', 'equal sd')                |                  0 |         262 |                8 |                  2 |              1 |                0 |
| ('senc', 2.0, 'national', 'interactions only')       |                  0 |         783 |             3097 |               4328 |           4359 |             5058 |
| ('redistrict', 0.5, 'indep', 'equal sd')             |                  0 |           0 |                0 |                  0 |              0 |                0 |
| ('redistrict', 0.5, 'national', 'equal sd')          |                  0 |           2 |                0 |                  0 |              0 |                0 |
| ('redistrict', 0.5, 'national', 'interactions only') |               1036 |           1 |             2927 |               4118 |           3320 |             1398 |
| ('redistrict', 2.0, 'indep', 'equal sd')             |                  0 |          36 |                1 |                  1 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'equal sd')          |                  0 |          40 |                2 |                  2 |              0 |                0 |
| ('redistrict', 2.0, 'national', 'interactions only') |                 10 |         363 |             2268 |               2926 |           2944 |              939 |

## Volume (middle area, sigma_b 1): median ESS per 1000 gradients

|                                   |   ILR + log volume |   log cells |   position logit |   odds-ratio logit |   adjusted row |   adjusted table |
|:----------------------------------|-------------------:|------------:|-----------------:|-------------------:|---------------:|-----------------:|
| ('redistrict', 5.23)              |             103.81 |       63.49 |            61.25 |              31.84 |          31.21 |            30.77 |
| ('redistrict', 523.0000000000001) |               7.16 |       39.13 |           145.59 |             130.39 |         143.21 |            42.27 |
| ('redistrict', 52300.0)           |               0.8  |        2.51 |           181.17 |             163.14 |         119.18 |            55.94 |
| ('senc', 15.9)                    |              54.71 |       43.86 |            56.42 |              33.61 |          41.74 |            30.03 |
| ('senc', 1590.0)                  |               6.12 |        3.94 |            87.54 |              76.22 |          64.33 |            47.76 |
| ('senc', 159000.0)                |               0.58 |        0.65 |            94.88 |              81.61 |          76.74 |            50.19 |

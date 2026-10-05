# Check that each allocation function works

12 random margins and parameter vectors per row (parameters ~ N(0, 1.5)), table total 1000; 200 random interior tables for the onto check. Orders: natural = rows in order, last row and last column by subtraction; shuffled = a random row order and a random subtraction column in each row.

| Scheme | Size | Order | Margins | Smallest cell | Inverse | Jacobian | Onto | Unreachable |
|---|---|---|---|---|---|---|---|---|
| between bounds, sharp | 2x2 | natural | 2.8e-16 | 1.1e-02 | 5.3e-15 | 1.3e-09 | 1.7e-16 | 0.0% |
| between bounds, sharp | 3x3 | natural | 2.3e-16 | 3.6e-04 | 2.7e-15 | 2.3e-09 | 1.7e-16 | 0.0% |
| between bounds, sharp | 3x3 | shuffled | 3.3e-16 | 5.6e-03 | 6.7e-15 | 2.8e-09 | 1.7e-16 | 0.0% |
| between bounds, sharp | 3x5 | natural | 3.4e-16 | 1.5e-04 | 1.6e-13 | 4.4e-09 | 1.7e-16 | 0.0% |
| between bounds, sharp | 3x5 | shuffled | 2.8e-16 | 3.8e-05 | 7.4e-13 | 6.4e-09 | 1.7e-16 | 0.0% |
| between bounds, sharp | 5x3 | natural | 2.4e-16 | 4.4e-05 | 1.4e-13 | 3.7e-09 | 1.4e-16 | 0.0% |
| between bounds, sharp | 5x3 | shuffled | 2.0e-16 | 3.3e-04 | 2.4e-14 | 5.6e-09 | 1.6e-16 | 0.0% |
| between bounds, sharp | 5x5 | natural | 2.6e-16 | 4.1e-05 | 1.4e-13 | 6.5e-09 | 1.9e-16 | 0.0% |
| between bounds, sharp | 5x5 | shuffled | 1.8e-16 | 3.5e-05 | 7.2e-14 | 5.7e-09 | 1.1e-16 | 0.0% |
| between bounds, sharp | 7x7 | natural | 3.7e-16 | 6.6e-07 | 3.7e-12 | 1.0e-08 | 1.3e-16 | 0.0% |
| between bounds, sharp | 7x7 | shuffled | 3.4e-16 | 8.9e-07 | 8.0e-13 | 9.4e-09 | 1.1e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 2x2 | natural | 2.3e-16 | 2.6e-03 | 8.5e-16 | 1.6e-09 | 1.7e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 3x3 | natural | 2.3e-16 | 1.4e-03 | 1.2e-14 | 2.0e-09 | 2.3e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 3x3 | shuffled | 2.8e-16 | 2.0e-04 | 2.8e-14 | 2.4e-09 | 3.4e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 3x5 | natural | 2.8e-16 | 6.4e-04 | 2.4e-14 | 4.2e-09 | 3.1e-16 | 0.5% |
| between bounds, smoothed (delta 0.05) | 3x5 | shuffled | 1.7e-16 | 1.4e-04 | 3.5e-14 | 2.8e-09 | 2.1e-16 | 0.5% |
| between bounds, smoothed (delta 0.05) | 5x3 | natural | 2.8e-16 | 4.0e-04 | 1.2e-14 | 3.8e-09 | 2.4e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 5x3 | shuffled | 1.7e-16 | 1.3e-04 | 1.8e-14 | 2.9e-09 | 2.4e-16 | 0.5% |
| between bounds, smoothed (delta 0.05) | 5x5 | natural | 1.4e-16 | 5.4e-06 | 2.7e-12 | 5.4e-09 | 1.6e-16 | 0.5% |
| between bounds, smoothed (delta 0.05) | 5x5 | shuffled | 2.8e-16 | 6.3e-06 | 1.3e-12 | 4.5e-09 | 2.0e-16 | 1.0% |
| between bounds, smoothed (delta 0.05) | 7x7 | natural | 3.4e-16 | 4.7e-06 | 5.0e-12 | 9.1e-09 | 1.8e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 7x7 | shuffled | 1.7e-16 | 8.0e-07 | 3.5e-13 | 4.9e-08 | 2.3e-16 | 0.5% |
| cell logit | 2x2 | natural | 1.7e-16 | 2.1e-03 | 3.9e-15 | 8.2e-10 | 8.0e-16 | 0.0% |
| cell logit | 3x3 | natural | 2.3e-16 | 2.8e-04 | 1.0e-13 | 1.9e-09 | 6.3e-16 | 0.0% |
| cell logit | 3x3 | shuffled | 2.8e-16 | 3.1e-03 | 6.7e-15 | 1.3e-09 | 1.7e-15 | 0.0% |
| cell logit | 3x5 | natural | 2.8e-16 | 1.6e-04 | 8.1e-14 | 3.5e-09 | 3.6e-16 | 0.0% |
| cell logit | 3x5 | shuffled | 2.3e-16 | 1.2e-04 | 2.3e-14 | 2.7e-09 | 4.8e-16 | 0.0% |
| cell logit | 5x3 | natural | 3.1e-16 | 1.9e-04 | 8.0e-15 | 2.8e-09 | 6.0e-16 | 0.0% |
| cell logit | 5x3 | shuffled | 2.8e-16 | 4.5e-05 | 5.1e-14 | 2.4e-09 | 1.2e-15 | 0.0% |
| cell logit | 5x5 | natural | 2.4e-16 | 4.1e-05 | 1.2e-12 | 3.8e-09 | 6.9e-15 | 0.0% |
| cell logit | 5x5 | shuffled | 2.3e-16 | 3.4e-06 | 3.8e-13 | 5.8e-09 | 4.8e-16 | 0.0% |
| cell logit | 7x7 | natural | 3.4e-16 | 4.5e-06 | 6.5e-14 | 7.0e-09 | 3.0e-15 | 0.0% |
| cell logit | 7x7 | shuffled | 2.7e-16 | 1.5e-05 | 4.2e-13 | 1.1e-08 | 9.9e-16 | 0.0% |
| adjusted row | 2x2 | natural | 2.3e-16 | 9.6e-03 | 8.9e-16 | 4.7e-10 | 2.6e-16 | 0.0% |
| adjusted row | 3x3 | natural | 2.8e-16 | 1.8e-03 | 1.4e-14 | 1.5e-09 | 1.6e-16 | 0.0% |
| adjusted row | 3x3 | shuffled | 1.7e-16 | 4.2e-04 | 9.3e-15 | 1.5e-09 | 1.6e-16 | 0.0% |
| adjusted row | 3x5 | natural | 1.1e-16 | 1.6e-04 | 2.7e-14 | 3.4e-09 | 9.9e-17 | 0.0% |
| adjusted row | 3x5 | shuffled | 2.8e-16 | 1.5e-04 | 3.5e-14 | 1.5e-09 | 1.1e-16 | 0.0% |
| adjusted row | 5x3 | natural | 3.7e-16 | 7.7e-05 | 1.5e-14 | 2.9e-09 | 2.0e-16 | 0.0% |
| adjusted row | 5x3 | shuffled | 2.9e-16 | 1.1e-04 | 3.0e-15 | 3.6e-09 | 1.7e-16 | 0.0% |
| adjusted row | 5x5 | natural | 2.3e-16 | 3.9e-06 | 3.6e-15 | 5.9e-09 | 1.7e-16 | 0.0% |
| adjusted row | 5x5 | shuffled | 2.8e-16 | 1.2e-04 | 3.3e-13 | 8.4e-09 | 1.2e-16 | 0.0% |
| adjusted row | 7x7 | natural | 3.4e-16 | 6.4e-08 | 2.1e-12 | 2.4e-08 | 1.1e-16 | 0.0% |
| adjusted row | 7x7 | shuffled | 3.7e-16 | 5.9e-05 | 1.1e-13 | 2.0e-08 | 7.1e-17 | 0.0% |
| adjusted table | 2x2 | natural | 3.4e-16 | 2.0e-02 | 1.3e-15 | 5.0e-10 | 1.0e-15 | 0.0% |
| adjusted table | 3x3 | natural | 4.5e-16 | 1.1e-04 | 1.6e-15 | 3.3e-09 | 9.1e-16 | 0.0% |
| adjusted table | 3x5 | natural | 5.7e-16 | 1.3e-04 | 1.4e-15 | 2.2e-09 | 5.1e-16 | 0.0% |
| adjusted table | 5x3 | natural | 4.5e-16 | 1.8e-04 | 1.8e-15 | 6.4e-09 | 6.3e-16 | 0.0% |
| adjusted table | 5x5 | natural | 5.1e-16 | 3.3e-05 | 1.6e-15 | 9.1e-09 | 2.6e-16 | 0.0% |
| adjusted table | 7x7 | natural | 2.4e-16 | 1.3e-06 | 1.3e-15 | 1.7e-08 | 1.8e-16 | 0.0% |

## Worst case by scheme

| Scheme | Margins | Smallest cell | Inverse | Jacobian | Onto | Unreachable |
|---|---|---|---|---|---|---|
| between bounds, sharp | 3.7e-16 | 6.6e-07 | 3.7e-12 | 1.0e-08 | 1.9e-16 | 0.0% |
| between bounds, smoothed (delta 0.05) | 3.4e-16 | 8.0e-07 | 5.0e-12 | 4.9e-08 | 3.4e-16 | 1.0% |
| cell logit | 3.4e-16 | 3.4e-06 | 1.2e-12 | 1.1e-08 | 6.9e-15 | 0.0% |
| adjusted row | 3.7e-16 | 6.4e-08 | 2.1e-12 | 2.4e-08 | 2.6e-16 | 0.0% |
| adjusted table | 5.7e-16 | 1.3e-06 | 1.8e-15 | 1.7e-08 | 1.0e-15 | 0.0% |

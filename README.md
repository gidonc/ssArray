# ssArray

Allocation functions for tables with fixed margins: smooth one-to-one maps from unconstrained
parameters onto the set of R x C tables with given row and column totals, for use as parameter
transforms in Stan. Research code and write-up; the R package that uses these ideas is
[ssEI](https://github.com/gidonc/ssEI).

## The four allocation functions

| File in `stan/functions/` | Name | How a table is built |
|---|---|---|
| `alloc_bounds.stan` | Between bounds | Cell by cell; each cell at a fraction of the way between its lower and upper bound |
| `alloc_cell_logit.stan` | Cell logit | Cell by cell; each cell placed by the log odds ratio of its collapsed 2 x 2 table |
| `alloc_row.stan` | Adjusted row | Row by row; each cell takes `capacity * inv_logit(eta + tau)`, one solved adjustment `tau` per row |
| `alloc_table.stan` | Adjusted table | Whole table at once, `exp(theta + u + v)`, row and column adjustments solved jointly (raking style) |

Each file has a forward function returning the table and the log Jacobian, and an inverse.

## Layout

- `stan/functions/` the allocation functions and the single-area parameter map (the single source of truth)
- `stan/single_area.stan` one area, one density on tables, sampled in a choice of six coordinate systems
- `stan/check_*.stan` evaluate the functions at supplied points
- `data/` example tables
- `experiments/` one script per result; outputs in `experiments/output/`
- `R/` R wrappers
- `paper/` the write-up

## Status

- `experiments/01_allocation_works.py`: all four functions meet the margins, invert, match a
  finite-difference Jacobian and reach every interior table, at 2x2 to 7x7.
- `experiments/02_single_area_models.py`: the single-area model with R x C parameters in six
  coordinate systems (ILR + log volume, log cells, and the four allocation functions with the margins
  as parameters). Jacobians checked; the sequential versions give the same posterior.
- The R files in `R/` have not been run.
- Experiments are in Python for now (run with cmdstanpy and CmdStan 2.32.2).

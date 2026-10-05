# Allocation functions for tables with fixed margins: outline

Scope of this first paper: one R x C table with both margins fixed. Third margins and
cross-area totals come later.

1. **The problem.** Tables with fixed margins as a parameter space of (R-1)(C-1) dimensions;
   the need for a smooth one-to-one map onto it with a tractable Jacobian.
2. **Sequential allocation works.** Any fill order with exact bounds reaches every feasible table
   exactly once; in the two-way case the bounds are a max and a min of two terms; the Jacobian is
   triangular. Checked numerically in `experiments/01_allocation_works.py`.
3. **Four allocation functions.** Between bounds, cell logit, adjusted row, adjusted table: each
   defined by its coordinates, its Jacobian and its cost.
4. **Curvature.** How far each function's coordinates are from linear in the within-row
   log-ratios, where the kinks are, and what that does to a sampler.
5. **What moves the curvature.** Table size, how extreme the table is, small cells, the choice
   of order and of the cell found by subtraction.
6. **Cost against geometry.**
7. **Bridge.** Structural zeros in a two-way table as the first case where the simple bounds are
   not exact and conditions have to be added.

## Naming

"Importance weights" was the working name for the between-bounds scheme. In Chen, Diaconis,
Holmes and Liu (2005) the importance weight is the correction for sampling a table from a
proposal; here nothing is sampled and the Jacobian plays that role. The write-up uses
"between bounds" for the scheme and cites that paper for the sequential construction.

# Make the Stan allocation functions callable from R.
#
# NOT YET RUN: written without access to R. Please check it and tell me what breaks.
# Needs cmdstanr (with CmdStan installed) and Rcpp.
#
# Usage, from the repository root:
#   source("R/expose.R")
#   ss_expose()
#   a <- alloc_row(w = c(300, 500, 200), m = c(400, 350, 250), lam = rep(0, 4),
#                  row_order = 1:3, rem_col = c(3L, 3L, 3L))
#   a[1:3, ]    # the table;  a[4, 1] is the log Jacobian
#
# Functions exposed (see the comments at the top of each file in stan/functions/):
#   alloc_bounds(w, m, lam, row_order, rem_col, delta), alloc_bounds_inv(T, row_order, rem_col, delta)
#   alloc_cell_logit(w, m, lam, row_order, rem_col),    alloc_cell_logit_inv(T, row_order, rem_col)
#   alloc_row(w, m, lam, row_order, rem_col),           alloc_row_inv(T, row_order, rem_col)
#   alloc_table(w, m, lam),                             alloc_table_inv(T)

ss_expose <- function(stan_dir = "stan") {
  mod <- cmdstanr::cmdstan_model(
    file.path(stan_dir, "check_alloc.stan"),
    include_paths = stan_dir,
    compile_standalone = TRUE
  )
  mod$expose_functions(global = TRUE)
  invisible(mod)
}

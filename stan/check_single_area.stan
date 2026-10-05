// Evaluates the single-area parameter-to-table map at supplied points (run with fixed_param).
functions {
#include functions/alloc_bounds.stan
#include functions/alloc_cell_logit.stan
#include functions/alloc_row.stan
#include functions/alloc_table.stan
#include functions/single_area.stan
}
data {
  int<lower=2> R;
  int<lower=2> C;
  int<lower=0, upper=5> param;
  vector<lower=0>[R] obs_w;
  vector<lower=0>[C] obs_m;
  real<lower=0> eps;
  int<lower=0, upper=1> scale_margins;
  matrix[R * C, R * C - 1] V;
  array[R] int<lower=1, upper=R> row_order;
  array[R] int<lower=1, upper=C> rem_col;
  real<lower=0> delta;
  int<lower=0> N;
  array[N] vector[R * C] theta;
}
generated quantities {
  array[N] matrix[R, C] T;
  vector[N] lj;
  for (n in 1:N) {
    matrix[R + 1, C] a = sa_table(theta[n], param, R, C, V, obs_w, obs_m, eps, scale_margins, row_order, rem_col, delta);
    T[n] = a[1:R, 1:C];
    lj[n] = a[R + 1, 1];
  }
}

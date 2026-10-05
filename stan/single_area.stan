// ONE AREA, one density on R x C tables, sampled in a choice of coordinates.
//
// Density on the table (the same whatever `param` is):
//   prior     ILR coordinates of the table ~ normal(mu_b, sigma_b), log total ~ normal(mu_logv, sigma_logv)
//             (a placeholder: substantive priors for real tables go here later)
//   margins   observed row and column totals ~ normal(the table's totals, eps * sqrt(observed total)),
//             so eps = 1 is about Poisson noise and smaller eps is a tighter penalty.
//
// `param` chooses the coordinates (see functions/single_area.stan):
//   0 = ILR + log volume (no allocation function), 5 = log cells,
//   1 = between bounds, 2 = cell logit, 3 = adjusted row, 4 = adjusted table (margins as parameters).
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
  vector<lower=0>[R] obs_w;                            // observed row totals
  vector<lower=0>[C] obs_m;                            // observed column totals
  real<lower=0> eps;                                   // relative width of the margin penalty
  matrix[R * C, R * C - 1] V;                          // ILR basis (orthonormal columns, each summing to zero), cells row by row
  vector[R * C - 1] mu_b;
  vector<lower=0>[R * C - 1] sigma_b;
  real mu_logv;
  real<lower=0> sigma_logv;
  array[R] int<lower=1, upper=R> row_order;            // sequential schemes: fill order of the rows
  array[R] int<lower=1, upper=C> rem_col;              // sequential schemes: column found by subtraction / reference
  real<lower=0> delta;                                 // between bounds: kink smoothing (0 = sharp)
}
parameters {
  vector[R * C] theta;
}
transformed parameters {
  matrix[R, C] T;
  real lj;
  {
    matrix[R + 1, C] a = sa_table(theta, param, R, C, V, obs_w, obs_m, eps, row_order, rem_col, delta);
    T = a[1:R, 1:C];
    lj = a[R + 1, 1];
  }
}
model {
  vector[R * C] lt = log(to_vector(T'));               // log cells, row by row
  // prior, stated on (ILR, log total); -sum(lt) takes it to a density on the cells, lj to a density on theta
  target += normal_lpdf(V' * lt | mu_b, sigma_b) + normal_lpdf(log(sum(T)) | mu_logv, sigma_logv) - sum(lt) + lj;
  // margin penalty
  for (r in 1:R) target += normal_lpdf(obs_w[r] | sum(T[r]), eps * sqrt(obs_w[r]));
  for (c in 1:C) target += normal_lpdf(obs_m[c] | sum(col(T, c)), eps * sqrt(obs_m[c]));
}
generated quantities {
  matrix[R, C] log_T = log(T);                         // for summaries: closer to normal than the cells
}

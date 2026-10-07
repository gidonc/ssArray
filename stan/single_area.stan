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
  int<lower=0, upper=1> scale_margins;                 // sequential versions: 1 = margin parameters on the scale of the penalty
  // The prior's basis. split_basis = 0: any orthonormal log-contrast basis, given as the matrix V (cells row by row);
  // the prior then costs a dense product of about (RC)^2 per gradient. split_basis = 1: the split basis (row effects,
  // column effects, interactions; see functions/single_area.stan), evaluated without its matrix at a cost of about
  // RC(R + C); V is then not used and is passed with zero rows and columns.
  int<lower=0, upper=1> split_basis;
  matrix[split_basis ? 0 : R * C, split_basis ? 0 : R * C - 1] V;
  vector[R * C - 1] mu_b;
  vector<lower=0>[R * C - 1] sigma_b;
  real mu_logv;
  real<lower=0> sigma_logv;
  array[R] int<lower=1, upper=R> row_order;            // sequential schemes: fill order of the rows
  array[R] int<lower=1, upper=C> rem_col;              // sequential schemes: column found by subtraction / reference
  real<lower=0> delta;                                 // between bounds: kink smoothing (0 = sharp)
  // sequential versions: the R + C - 1 margin parameters used by sa_table are K_margin * theta[1:(R + C - 1)].
  // A fixed linear change of margin coordinates (constant Jacobian, so the density on tables is unchanged):
  // identity = column log-ratios against the last column; other choices: against the largest column, or whitened
  // by the margin penalty so that no column is a reference (see sa_margin_K in R/single_area.R).
  matrix[R + C - 1, R + C - 1] K_margin;
  // adjusted row and adjusted table: 1 = the interior parameters are orthonormal contrasts with no reference cell.
  //   adjusted row    the C - 1 parameters of a row are Helmert contrasts among all C of its logits
  //                   (0 = logits against the row's reference column rem_col)
  //   adjusted table  the parameters are the interaction coordinates of the split basis
  //                   (0 = log odds ratios against the last row and last column)
  // Both are fixed linear changes of coordinates, so the density on tables is unchanged.
  int<lower=0, upper=1> centred_interior;
}
transformed data {
  matrix[R, R - 1] Hr = sa_helmert(R);
  matrix[C, C - 1] Hc = sa_helmert(C);
}
parameters {
  vector[R * C] theta;
}
transformed parameters {
  matrix[R, C] T;
  real lj;
  {
    int D = R * C;
    int nm = R + C - 1;
    vector[D] th = theta;
    if (param >= 1 && param <= 4) th[1:nm] = K_margin * theta[1:nm];
    if (centred_interior == 1 && param == 3) {
      for (k in 1:(R - 1)) {
        int last = rem_col[row_order[k]];
        int off = nm + (k - 1) * (C - 1);
        vector[C] e = Hc * theta[(off + 1):(off + C - 1)];
        int i = 0;
        for (c in 1:C) if (c != last) { i += 1; th[off + i] = e[c] - e[last]; }
      }
    }
    if (centred_interior == 1 && param == 4) {
      matrix[R, C] I = Hr * to_matrix(theta[(nm + 1):D], C - 1, R - 1)' * Hc';
      for (r in 1:(R - 1)) for (c in 1:(C - 1))
        th[nm + (r - 1) * (C - 1) + c] = I[r, c] - I[r, C] - I[R, c] + I[R, C];
    }
    if (param == 0 && split_basis == 1) {
      matrix[R, C] L = split_logcells(theta[1:(D - 1)], Hr, Hc);
      L += theta[D] - log_sum_exp(L);
      T = exp(L);
      lj = sum(L);                                     // as in sa_table: log cells are linear in (coordinates, mean log cell)
    } else {
      matrix[R + 1, C] a = sa_table(th, param, R, C, V, obs_w, obs_m, eps, scale_margins, row_order, rem_col, delta);
      T = a[1:R, 1:C];
      lj = a[R + 1, 1];
    }
  }
}
model {
  vector[R * C] lt = log(to_vector(T'));               // log cells, row by row
  // prior, stated on (ILR, log total); -sum(lt) takes it to a density on the cells, lj to a density on theta
  if (split_basis == 1) target += normal_lpdf(split_coords(log(T), Hr, Hc) | mu_b, sigma_b);
  else target += normal_lpdf(V' * lt | mu_b, sigma_b);
  target += normal_lpdf(log(sum(T)) | mu_logv, sigma_logv) - sum(lt) + lj;
  // margin penalty
  for (r in 1:R) target += normal_lpdf(obs_w[r] | sum(T[r]), eps * sqrt(obs_w[r]));
  for (c in 1:C) target += normal_lpdf(obs_m[c] | sum(col(T, c)), eps * sqrt(obs_m[c]));
}
generated quantities {
  matrix[R, C] log_T = log(T);                         // for summaries: closer to normal than the cells
}

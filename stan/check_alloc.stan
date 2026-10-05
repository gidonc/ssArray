// Evaluates one allocation function at supplied points (no sampling: run with fixed_param).
//   scheme 1 = between bounds, 2 = cell logit, 3 = adjusted row, 4 = adjusted table
functions {
#include functions/alloc_bounds.stan
#include functions/alloc_cell_logit.stan
#include functions/alloc_row.stan
#include functions/alloc_table.stan
}
data {
  int<lower=2> R;
  int<lower=2> C;
  int<lower=1, upper=4> scheme;
  array[R] int<lower=1, upper=R> row_order;
  array[R] int<lower=1, upper=C> rem_col;
  real<lower=0> delta;                                 // kink smoothing, scheme 1 only
  int<lower=0> N;                                      // forward evaluations
  array[N] vector[R] w;
  array[N] vector[C] m;
  array[N] vector[(R - 1) * (C - 1)] lam;
  int<lower=0> N2;                                     // tables to invert and rebuild
  array[N2] matrix[R, C] Tin;
}
generated quantities {
  array[N] matrix[R, C] T;
  vector[N] lj;
  array[N] vector[(R - 1) * (C - 1)] lam_back;
  array[N2] vector[(R - 1) * (C - 1)] lam_in;
  array[N2] matrix[R, C] T_back;
  for (n in 1:N) {
    matrix[R + 1, C] a;
    if (scheme == 1) a = alloc_bounds(w[n], m[n], lam[n], row_order, rem_col, delta);
    else if (scheme == 2) a = alloc_cell_logit(w[n], m[n], lam[n], row_order, rem_col);
    else if (scheme == 3) a = alloc_row(w[n], m[n], lam[n], row_order, rem_col);
    else a = alloc_table(w[n], m[n], lam[n]);
    T[n] = a[1:R, 1:C];
    lj[n] = a[R + 1, 1];
    if (scheme == 1) lam_back[n] = alloc_bounds_inv(T[n], row_order, rem_col, delta);
    else if (scheme == 2) lam_back[n] = alloc_cell_logit_inv(T[n], row_order, rem_col);
    else if (scheme == 3) lam_back[n] = alloc_row_inv(T[n], row_order, rem_col);
    else lam_back[n] = alloc_table_inv(T[n]);
  }
  for (n in 1:N2) {
    vector[R] wn;
    vector[C] mn;
    matrix[R + 1, C] a;
    for (r in 1:R) wn[r] = sum(Tin[n][r]);
    for (c in 1:C) mn[c] = sum(col(Tin[n], c));
    if (scheme == 1) lam_in[n] = alloc_bounds_inv(Tin[n], row_order, rem_col, delta);
    else if (scheme == 2) lam_in[n] = alloc_cell_logit_inv(Tin[n], row_order, rem_col);
    else if (scheme == 3) lam_in[n] = alloc_row_inv(Tin[n], row_order, rem_col);
    else lam_in[n] = alloc_table_inv(Tin[n]);
    if (scheme == 1) a = alloc_bounds(wn, mn, lam_in[n], row_order, rem_col, delta);
    else if (scheme == 2) a = alloc_cell_logit(wn, mn, lam_in[n], row_order, rem_col);
    else if (scheme == 3) a = alloc_row(wn, mn, lam_in[n], row_order, rem_col);
    else a = alloc_table(wn, mn, lam_in[n]);
    T_back[n] = a[1:R, 1:C];
  }
}

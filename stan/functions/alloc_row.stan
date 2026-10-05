// ---------------------------------------------------------------------------------------------
// ADJUSTED-ROW allocation of an R x C table with fixed margins.
// Rows are filled in the order row_order[1:(R-1)]; row_order[R] is found by subtraction.
// A row is allocated across ALL its columns in one step:
//
//     x[c] = cap[c] * inv_logit(eta[c] + tau),
//
//   cap[c] = what earlier rows have left in column c,
//   eta[c] = the row's free logits (eta[rem_col[r]] = 0, the reference column),
//   tau    = the one adjustment per row: the common shift that makes the row add up to its total.
//
// The row total is increasing in tau, so tau is a one-dimensional monotone solve.
// Every cell is strictly between 0 and its column's capacity, the row total is met exactly, the whole
// feasible range is reachable and there are no kinks.
// ---------------------------------------------------------------------------------------------

// the shift tau with  sum_c cap[c] * inv_logit(eta[c] + tau) = total   (safeguarded Newton)
real row_shift(vector cap, vector eta, real total) {
  real base = logit(total / sum(cap));
  real lo = base - max(eta);                       // row total at lo <= total
  real hi = base - min(eta);                       // row total at hi >= total
  real tau = fmin(fmax(log(total) - log_sum_exp(log(cap) + eta), lo), hi);   // exact when the fractions are small
  for (it in 1:80) {
    vector[rows(cap)] p = inv_logit(eta + tau);
    real f = dot_product(cap, p) - total;
    real fp = dot_product(cap, p .* (1 - p));
    real nt = tau - f / fp;
    // once converged, one more Newton step is exact to rounding error and carries the right gradient
    if (abs(f) < 1e-10 * total) { tau = nt; break; }
    if (f > 0) hi = tau; else lo = tau;
    tau = (nt > lo && nt < hi) ? nt : 0.5 * (lo + hi);
  }
  return tau;
}

// lam: (R-1)*(C-1) logits, row by row in allocation order, C-1 per row (columns other than rem_col[r]).
// Returns (R+1) x C: rows 1:R = the table, [R+1, 1] = log |d free cells / d lam|.
matrix alloc_row(vector w, vector m, vector lam, array[] int row_order, array[] int rem_col) {
  int R = rows(w);
  int C = rows(m);
  vector[C] sc = m;
  matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
  real lj = 0;
  for (k in 1:(R - 1)) {
    int r = row_order[k];
    int last = rem_col[r];
    int i = 0;
    vector[C] eta = rep_vector(0, C);
    vector[C] x;
    vector[C] v;
    real tau;
    for (c in 1:C) if (c != last) { i += 1; eta[c] = lam[(k - 1) * (C - 1) + i]; }
    tau = row_shift(sc, eta, w[r]);
    x = sc .* inv_logit(eta + tau);
    v = x .* (1 - x ./ sc);                          // d x[c] / d (eta[c] + tau)
    for (c in 1:C) out[r, c] = x[c];
    lj += sum(log(v)) - log(sum(v));                 // rank-one update: the row's Jacobian in closed form
    sc -= x;
  }
  for (c in 1:C) out[row_order[R], c] = sc[c];
  out[R + 1, 1] = lj;
  return out;
}

vector alloc_row_inv(matrix T, array[] int row_order, array[] int rem_col) {
  int R = rows(T);
  int C = cols(T);
  vector[C] sc;
  vector[(R - 1) * (C - 1)] lam;
  for (c in 1:C) sc[c] = sum(col(T, c));
  for (k in 1:(R - 1)) {
    int r = row_order[k];
    int last = rem_col[r];
    int i = 0;
    for (c in 1:C) if (c != last) {
      i += 1;
      lam[(k - 1) * (C - 1) + i] = logit(T[r, c] / sc[c]) - logit(T[r, last] / sc[last]);
    }
    sc -= T[r]';
  }
  return lam;
}

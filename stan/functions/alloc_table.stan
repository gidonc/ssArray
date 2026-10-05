// ---------------------------------------------------------------------------------------------
// ADJUSTED-TABLE allocation of an R x C table with fixed margins (raking style).
// The whole table is allocated at once:
//
//     T[r, c] = exp(theta[r, c] + u[r] + v[c]),
//
//   theta = the free parameters: (R-1) x (C-1) log odds ratios against the last row and last column
//           (theta is 0 in the last row and the last column),
//   u, v  = one adjustment per row and per column, solved jointly so that both margins are met.
//
// This is the table that raking (iterative proportional fitting) of the seed exp(theta) converges to.
// u and v minimise a convex function, so the solution exists and is unique for any theta; it is found by
// a few raking sweeps followed by Newton steps. There is no fill order and there are no kinks.
// ---------------------------------------------------------------------------------------------

// lam: (R-1)*(C-1) log odds ratios, row by row (rows 1:(R-1), columns 1:(C-1)).
// Returns (R+1) x C: rows 1:R = the table, [R+1, 1] = log |d free cells / d lam|.
matrix alloc_table(vector w, vector m, vector lam) {
  int R = rows(w);
  int C = rows(m);
  int K = (R - 1) * (C - 1);
  int D = R + C - 1;
  real tot = sum(w);
  matrix[R, C] th = rep_matrix(0, R, C);
  vector[R] u;
  vector[C] v = log(m) - log(tot);
  matrix[R, C] T;
  matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
  for (r in 1:(R - 1)) for (c in 1:(C - 1)) th[r, c] = lam[(r - 1) * (C - 1) + c];
  // raking sweeps to get close
  for (s in 1:3) {
    for (r in 1:R) u[r] = log(w[r]) - log_sum_exp(th[r]' + v);
    for (c in 1:C) v[c] = log(m[c]) - log_sum_exp(col(th, c) + u);
  }
  u += v[C];
  v -= v[C];                                         // v[C] = 0 fixes the one redundant direction
  // Newton on (u, v[1:(C-1)])
  for (it in 1:50) {
    vector[D] g;
    matrix[D, D] H = rep_matrix(0, D, D);
    vector[D] step;
    real big;
    for (r in 1:R) for (c in 1:C) T[r, c] = exp(th[r, c] + u[r] + v[c]);
    for (r in 1:R) { g[r] = sum(T[r]) - w[r]; H[r, r] = sum(T[r]); }
    for (c in 1:(C - 1)) {
      g[R + c] = sum(col(T, c)) - m[c];
      H[R + c, R + c] = sum(col(T, c));
      for (r in 1:R) { H[r, R + c] = T[r, c]; H[R + c, r] = T[r, c]; }
    }
    step = -(H \ g);
    big = max(abs(step));
    if (big > 2) step *= 2 / big;                    // safeguard far from the solution
    u += step[1:R];
    for (c in 1:(C - 1)) v[c] += step[R + c];
    // once converged, the step just taken is exact to rounding error and carries the right gradient
    if (max(abs(g)) < 1e-11 * tot) break;
  }
  for (r in 1:R) for (c in 1:C) T[r, c] = exp(th[r, c] + u[r] + v[c]);
  // Jacobian. The inverse map is explicit, theta[r, c] = log T[r, c] + log T[R, C] - log T[r, C] - log T[R, c],
  // with the last row and column linear in the free cells, so
  //   d theta[r, c] / d T[r', c'] = (r = r')(c = c') / T[r, c] + 1 / T[R, C] + (r = r') / T[r, C] + (c = c') / T[R, c].
  {
    matrix[K, K] J;
    for (r in 1:(R - 1)) for (c in 1:(C - 1)) {
      int a = (r - 1) * (C - 1) + c;
      for (r2 in 1:(R - 1)) for (c2 in 1:(C - 1)) {
        int b = (r2 - 1) * (C - 1) + c2;
        J[a, b] = inv(T[R, C]) + (r == r2 ? inv(T[r, C]) : 0) + (c == c2 ? inv(T[R, c]) : 0) + (a == b ? inv(T[r, c]) : 0);
      }
    }
    out[R + 1, 1] = -log_determinant_spd(J);
  }
  out[1:R, 1:C] = T;
  return out;
}

vector alloc_table_inv(matrix T) {
  int R = rows(T);
  int C = cols(T);
  vector[(R - 1) * (C - 1)] lam;
  for (r in 1:(R - 1)) for (c in 1:(C - 1))
    lam[(r - 1) * (C - 1) + c] = log(T[r, c]) + log(T[R, C]) - log(T[r, C]) - log(T[R, c]);
  return lam;
}

// ---------------------------------------------------------------------------------------------
// ONE AREA, R x C parameters: the map from an unconstrained parameter vector to the table.
// Every parameterisation has R * C parameters and maps one-to-one onto the positive R x C tables,
// so the same density on tables can be sampled in any of them.
//
//   param 0 : ILR + log volume.   theta[1:(RC-1)] = ILR coordinates in the basis V, theta[RC] = log total.
//   param 5 : log cell values, row by row.
//   param 1-4 : SEQUENTIAL.  theta = (z_row[R], z_col[C-1], lam[(R-1)(C-1)]):
//       row margins      w[r] = obs_w[r] * exp(s_r * z_row[r])
//       column margins   m    = sum(w) * softmax(log-ratios of obs_m + s_c .* z_col, 0)
//       scale_margins = 0 : s_r = s_c = 1, so z_row and z_col are plain log deviations from the observed margins
//       scale_margins = 1 : s_r = eps / sqrt(obs_w[r]), s_c = eps * sqrt(1/obs_m[c] + 1/obs_m[C]): the scale of a
//                           margin penalty of relative width eps, so z_row and z_col are roughly unit scale
//       interior         the table with margins (w, m) and allocation parameters lam, built by
//                        1 = between bounds, 2 = cell logit, 3 = adjusted row, 4 = adjusted table.
//     Either way the margins are parameters in their own right, so a margin penalty acts on its own coordinates.
//
// Needs the four alloc_*.stan files to be included first.
// Returns (R+1) x C: rows 1:R = the table, [R+1, 1] = log |d cells / d theta|.
// ---------------------------------------------------------------------------------------------
matrix sa_table(vector theta, int param, int R, int C, matrix V, vector obs_w, vector obs_m, real eps,
                int scale_margins, array[] int row_order, array[] int rem_col, real delta) {
  int D = R * C;
  matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
  if (param == 0) {
    vector[D] lt = V * theta[1:(D - 1)];
    lt += theta[D] - log_sum_exp(lt);
    for (r in 1:R) for (c in 1:C) out[r, c] = exp(lt[(r - 1) * C + c]);
    out[R + 1, 1] = sum(lt);                         // up to a constant: log cells are linear in (ILR, mean log cell)
  } else if (param == 5) {
    for (r in 1:R) for (c in 1:C) out[r, c] = exp(theta[(r - 1) * C + c]);
    out[R + 1, 1] = sum(theta);
  } else {
    vector[R] s_r = rep_vector(1, R);
    vector[C - 1] s_c = rep_vector(1, C - 1);
    if (scale_margins == 1) {
      s_r = eps * inv_sqrt(obs_w);
      s_c = eps * sqrt(inv(obs_m[1:(C - 1)]) + inv(obs_m[C]));
    }
    vector[R] w = obs_w .* exp(s_r .* theta[1:R]);
    real N = sum(w);
    vector[C] m = N * softmax(append_row(log(obs_m[1:(C - 1)]) - log(obs_m[C]) + s_c .* theta[(R + 1):(R + C - 1)], 0));
    vector[(R - 1) * (C - 1)] lam = theta[(R + C):D];
    if (param == 1) out = alloc_bounds(w, m, lam, row_order, rem_col, delta);
    else if (param == 2) out = alloc_cell_logit(w, m, lam, row_order, rem_col);
    else if (param == 3) out = alloc_row(w, m, lam, row_order, rem_col);
    else out = alloc_table(w, m, lam);
    // block triangular: (z_row -> w), (z_col -> m given the total), (lam -> interior given the margins)
    out[R + 1, 1] += sum(log(w)) + sum(log(s_r)) + sum(log(m)) - log(N) + sum(log(s_c));
  }
  return out;
}

// ---------------------------------------------------------------------------------------------
// The SPLIT basis without its matrix.
// Orthonormal log-contrast basis of an R x C table in three blocks: R-1 row effects, C-1 column effects and
// (R-1)(C-1) interactions, each built from Helmert contrasts (Hr, Hc). As a matrix it is
//   [ Hr (x) 1/sqrt(C),  1/sqrt(R) (x) Hc,  Hr (x) Hc ]        (cells row by row),
// which has about (RC)^2 entries. Because every block is a Kronecker product the coordinates can be had
// from the R x C matrix of log cells with small matrix products, at a cost of about RC(R + C).
// ---------------------------------------------------------------------------------------------
matrix sa_helmert(int n) {
  matrix[n, n - 1] H = rep_matrix(0, n, n - 1);
  for (k in 1:(n - 1)) {
    for (i in 1:k) H[i, k] = inv_sqrt(k * (k + 1.0));
    H[k + 1, k] = -k * inv_sqrt(k * (k + 1.0));
  }
  return H;
}

// coordinates of the log cells L: row effects, column effects, interactions (row by row)
vector split_coords(matrix L, matrix Hr, matrix Hc) {
  int R = rows(L);
  int C = cols(L);
  vector[R * C - 1] z;
  z[1:(R - 1)] = Hr' * (L * rep_vector(inv_sqrt(C), C));
  z[R:(R + C - 2)] = (rep_row_vector(inv_sqrt(R), R) * L * Hc)';
  z[(R + C - 1):(R * C - 1)] = to_vector((Hr' * L * Hc)');
  return z;
}

// the centred log cells with coordinates z (inverse of split_coords on tables whose log cells sum to zero)
matrix split_logcells(vector z, matrix Hr, matrix Hc) {
  int R = rows(Hr);
  int C = rows(Hc);
  return rep_matrix(Hr * z[1:(R - 1)] * inv_sqrt(C), C)
         + rep_matrix((Hc * z[R:(R + C - 2)])' * inv_sqrt(R), R)
         + Hr * to_matrix(z[(R + C - 1):(R * C - 1)], C - 1, R - 1)' * Hc';
}

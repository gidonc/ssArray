// ---------------------------------------------------------------------------------------------
// ONE AREA, R x C parameters: the map from an unconstrained parameter vector to the table.
// Every parameterisation has R * C parameters and maps one-to-one onto the positive R x C tables,
// so the same density on tables can be sampled in any of them.
//
//   param 0 : ILR + log volume.   theta[1:(RC-1)] = ILR coordinates in the basis V, theta[RC] = log total.
//   param 5 : log cell values, row by row.
//   param 1-4 : SEQUENTIAL.  theta = (z_row[R], z_col[C-1], lam[(R-1)(C-1)]):
//       row margins      w[r] = obs_w[r] * exp(s_r * z_row[r]),           s_r = eps / sqrt(obs_w[r])
//       column margins   m    = sum(w) * softmax(log-ratios of obs_m + s_c .* z_col, 0),
//                                                                  s_c = eps * sqrt(1/obs_m[c] + 1/obs_m[C])
//       interior         the table with margins (w, m) and allocation parameters lam, built by
//                        1 = between bounds, 2 = cell logit, 3 = adjusted row, 4 = adjusted table.
//     The margins are parameters in their own right, on the scale of a margin penalty of relative
//     width eps, so z_row and z_col are roughly unit scale however tight the penalty is.
//
// Needs the four alloc_*.stan files to be included first.
// Returns (R+1) x C: rows 1:R = the table, [R+1, 1] = log |d cells / d theta|.
// ---------------------------------------------------------------------------------------------
matrix sa_table(vector theta, int param, int R, int C, matrix V, vector obs_w, vector obs_m, real eps,
                array[] int row_order, array[] int rem_col, real delta) {
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
    vector[R] s_r = eps * inv_sqrt(obs_w);
    vector[C - 1] s_c = eps * sqrt(inv(obs_m[1:(C - 1)]) + inv(obs_m[C]));
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

// ---------------------------------------------------------------------------------------------
// BETWEEN-BOUNDS allocation of an R x C table with fixed margins.
// (The original ssEI scheme; "importance weights" in the working notes. The construction follows the
//  sequential fill of Chen, Diaconis, Holmes and Liu (2005); here each cell is placed by a parameter,
//  not sampled, and the Jacobian takes the place of their importance weight.)
//
// Rows are filled in the order row_order[1:(R-1)]; row_order[R] is found by subtraction.
// Within row r the columns are visited in natural order, skipping rem_col[r], which is found by subtraction.
// A visited cell has
//     lower = max(0, sr - rest),   upper = min(sc, sr),
//   sr   = what is left of the row,
//   sc   = what is left of this column,
//   rest = what is left of the columns this row has not visited yet (they must absorb the rest of the row),
// and is placed at  lower + inv_logit(lam) * (upper - lower).
//
// delta = 0 : sharp bounds (exact, kinked where sr = rest or sc = sr).
// delta > 0 : each kink is rounded over about delta x (the cell's range), inward only, so every table is
//             still feasible and the margins are still exact, but at most about 1.4 x delta of a cell's
//             range near a tie cannot be reached.
// ---------------------------------------------------------------------------------------------

// [lower, width] of one cell
vector bounds_cell(real sr, real sc, real rest, real delta) {
  vector[2] out;
  if (delta > 0) {
    real k = 1 / delta;
    real l_sr = log(fmax(sr, 1e-300));
    real l_sc = log(fmax(sc, 1e-300));
    real l_rest = log(fmax(rest, 1e-300));
    real l_oth = log(fmax(sc + rest - sr, 1e-300));
    // smooth version, from below, of the true range min(sc, sr, rest, sc + rest - sr)
    real l_w = -log_sum_exp([-k * l_sc, -k * l_sr, -k * l_rest, -k * l_oth]') / k;
    real ke = k * exp(fmin(l_sr - l_w, 30));
    real up = sc * exp(-log1p_exp(ke * (l_sc - l_sr)) / ke);
    real lo = sr * (-expm1(-log1p_exp(ke * (l_sr - l_rest)) / ke));
    out[1] = lo;
    out[2] = up - lo;
  } else {
    real lo = fmax(0, sr - rest);
    out[1] = lo;
    out[2] = fmin(sc, sr) - lo;
  }
  return out;
}

// lam: (R-1)*(C-1) logits, laid out row by row in allocation order, C-1 per row
// (the columns other than rem_col[r], in natural order).
// Returns (R+1) x C: rows 1:R = the table, [R+1, 1] = log |d free cells / d lam|.
matrix alloc_bounds(vector w, vector m, vector lam, array[] int row_order, array[] int rem_col, real delta) {
  int R = rows(w);
  int C = rows(m);
  vector[R] sr = w;
  vector[C] sc = m;
  matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
  real lj = 0;
  for (k in 1:(R - 1)) {
    int r = row_order[k];
    int last = rem_col[r];
    int i = 0;
    real rest = sum(sc);
    for (c in 1:C) if (c != last) {
      vector[2] b;
      real l;
      real x;
      i += 1;
      rest -= sc[c];
      b = bounds_cell(sr[r], sc[c], rest, delta);
      l = lam[(k - 1) * (C - 1) + i];
      x = b[1] + inv_logit(l) * b[2];
      out[r, c] = x;
      lj += log(b[2]) + log_inv_logit(l) + log1m_inv_logit(l);
      sc[c] -= x;
      sr[r] -= x;
    }
    out[r, last] = sr[r];
    sc[last] -= sr[r];
  }
  for (c in 1:C) out[row_order[R], c] = sc[c];
  out[R + 1, 1] = lj;
  return out;
}

// the lam that reproduce table T (on T's own margins)
vector alloc_bounds_inv(matrix T, array[] int row_order, array[] int rem_col, real delta) {
  int R = rows(T);
  int C = cols(T);
  vector[R] sr;
  vector[C] sc;
  vector[(R - 1) * (C - 1)] lam;
  for (r in 1:R) sr[r] = sum(T[r]);
  for (c in 1:C) sc[c] = sum(col(T, c));
  for (k in 1:(R - 1)) {
    int r = row_order[k];
    int last = rem_col[r];
    int i = 0;
    real rest = sum(sc);
    for (c in 1:C) if (c != last) {
      vector[2] b;
      i += 1;
      rest -= sc[c];
      b = bounds_cell(sr[r], sc[c], rest, delta);
      lam[(k - 1) * (C - 1) + i] = logit((T[r, c] - b[1]) / b[2]);
      sc[c] -= T[r, c];
      sr[r] -= T[r, c];
    }
    sc[last] -= sr[r];
  }
  return lam;
}

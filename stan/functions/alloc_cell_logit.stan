// ---------------------------------------------------------------------------------------------
// CELL-LOGIT allocation of an R x C table with fixed margins.
// Same sequential fill and order as alloc_bounds, but each visited cell x is placed by the log odds
// ratio of the 2 x 2 table it sits in once the rest of the row and the rest of the column are collapsed:
//
//        [ x        sr - x          ]   this row
//        [ sc - x   rest - sr + x   ]   the rows not yet filled
//          this col   columns not yet visited
//
//     lam = log x + log(rest - sr + x) - log(sr - x) - log(sc - x).
//
// x is the root of a quadratic. It lies strictly inside max(0, sr - rest) < x < min(sc, sr), covers that
// whole range as lam runs over the real line, and is smooth in sr, sc and rest: there are no kinks.
// lam = 0 is independence, x = sr * sc / (sc + rest).
// ---------------------------------------------------------------------------------------------

real cell_logit_cell(real sr, real sc, real rest, real lam) {
  real psi = exp(lam);
  real A = psi - 1;
  real B = -(psi * (sr + sc) + rest - sr);
  real Cq = psi * sr * sc;
  real disc = sqrt(fmax(B * B - 4 * A * Cq, 0));
  if (B <= 0) return 2 * Cq / (-B + disc);       // stable form of the root inside the bounds
  return (-B - disc) / (2 * A);
}

// log |dx / dlam| for that cell
real cell_logit_logjac(real x, real sr, real sc, real rest) {
  return -log(inv(x) + inv(sr - x) + inv(sc - x) + inv(rest - sr + x));
}

// Layout of lam and of the result as in alloc_bounds.
matrix alloc_cell_logit(vector w, vector m, vector lam, array[] int row_order, array[] int rem_col) {
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
      real x;
      i += 1;
      rest -= sc[c];
      x = cell_logit_cell(sr[r], sc[c], rest, lam[(k - 1) * (C - 1) + i]);
      out[r, c] = x;
      lj += cell_logit_logjac(x, sr[r], sc[c], rest);
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

vector alloc_cell_logit_inv(matrix T, array[] int row_order, array[] int rem_col) {
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
      real x = T[r, c];
      i += 1;
      rest -= sc[c];
      lam[(k - 1) * (C - 1) + i] = log(x) + log(rest - sr[r] + x) - log(sr[r] - x) - log(sc[c] - x);
      sc[c] -= x;
      sr[r] -= x;
    }
    sc[last] -= sr[r];
  }
  return lam;
}

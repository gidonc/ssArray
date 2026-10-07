"""Margin coordinates for the sequential versions of single_area.stan.

The model's margin parameters are z = (z_row[R], z_col[C-1]):  row totals obs_w * exp(s_r z_row), column shares a softmax of
(observed log-ratios against the LAST column + s_c z_col).  The data matrix K_margin sets z = K_margin * u, so u are the
sampled parameters.  Any K is a fixed linear change of coordinates: the density on tables is the same.
  'last'      K = identity (the original coordinates)
  'largest'   column log-ratios against the largest column instead of the last
  'whitened'  u has roughly unit, uncorrelated scale under the margin penalty (no reference column); use scale_margins = 0
margin_K returns (K, scale_margins to pass to the model).  Default 'largest' (as in R/single_area.R); the model's own default data
(base_data in 02) is still the identity, so earlier experiments reproduce.
"""
import numpy as np


def margin_K(obs_w, obs_m, eps, kind="largest", scale_margins=1):
    w = np.asarray(obs_w, float); m = np.asarray(obs_m, float); R, C = len(w), len(m); n = R + C - 1
    if kind == "last":
        return np.eye(n), scale_margins
    if kind == "largest":
        j = int(np.argmax(m))
        if j == C - 1: return np.eye(n), scale_margins
        s_old = eps * np.sqrt(1 / m[:C - 1] + 1 / m[C - 1]) if scale_margins else np.ones(C - 1)      # scale of the old column parameters
        others = [c for c in range(C) if c != j]                                                       # new parameters: columns other than j, natural order
        t_new = eps * np.sqrt(1 / m[others] + 1 / m[j]) if scale_margins else np.ones(C - 1)
        # e_c = log-ratio deviation of column c against column j; old d_c (against the last column) = e_c - e_last, d_j = -e_last
        A = np.zeros((C - 1, C - 1)); last = others.index(C - 1)
        for c in range(C - 1):
            if c != j: A[c, others.index(c)] += 1.0
            A[c, last] -= 1.0
        K = np.eye(n); K[R:, R:] = (A * t_new[None, :]) / s_old[:, None]
        return K, scale_margins
    if kind == "whitened":
        N = w.sum(); p = w / N; q = m / N
        M = np.zeros((R + C, n))                                    # d (log w, log m) / d z  with unit scales
        M[:R, :R] = np.eye(R)
        M[R:, :R] = p[None, :]                                      # the total moves every column
        B = np.zeros((C, C - 1)); B[:C - 1] = np.eye(C - 1); B = B - q[None, :C - 1]      # softmax: b_c - sum_c q_c b_c, b_last = 0
        M[R:, R:] = B
        H = M.T @ np.diag(np.r_[w, m] / eps ** 2) @ M + np.eye(n)   # penalty precision, + 1 so no parameter is wider than sd 1
        L = np.linalg.cholesky(H)
        return np.linalg.inv(L).T, 0
    raise ValueError(kind)

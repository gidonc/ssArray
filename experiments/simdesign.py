"""Simulated single-area designs with separate dials.  A design is what the model sees: margins, prior centre, prior sds.

SIZE
  R, C            table shape
  volume          (unit, value): which count is held fixed.  ('total', N), ('margin', mean count per row margin),
                  ('cell', mean count per cell).  They are one scale, N = margin * R = cell * R * C; the dial says which
                  one stays put when R and C change.
  eps             margin noise relative to count noise (sd = eps * sqrt(margin))
MARGINS (rows and columns separately)
  balance_row, balance_col   effective number of categories / number of categories, in (0, 1].  1 = equal margins.
  family_row, family_col     how imbalance is arranged: 'geometric' (spread evenly on a log scale),
                             'dominant' (one large margin, the rest equal), 'tiny' (one small margin, the rest equal)
  align           how the column ordering relates to the rows: 'same' (largest row meets largest column on the
                  diagonal), 'reversed', 'shifted' (columns rotated by one)
  jitter          log-normal sd applied to every margin, for replicates of the same kind
INTERIOR / PRIOR
  kappa           log odds ratio on the diagonal of the prior centre; the centre is raked to the margins
  sigma_b         prior sd on the interaction coordinates
  sigma_main      prior sd on the row and column effects (None = same as sigma_b; large = interactions-only prior)
  conflict        centre's row and column effects moved this many prior sds from the margins
make_design returns the data pieces and a dict of derived features, so a run can be read against what the design
actually is (counts per margin and per cell, effective numbers, tightness index, share of mass on the diagonal ...).
"""
import numpy as np


def helmert(D):
    V = np.zeros((D, D - 1))
    for k in range(1, D):
        V[:k, k - 1] = 1 / np.sqrt(k * (k + 1.0)); V[k, k - 1] = -k / np.sqrt(k * (k + 1.0))
    return V


def split_basis(R, C):
    """orthonormal log-contrast basis: R-1 row effects, C-1 column effects, then (R-1)(C-1) interactions"""
    one = lambda n: np.ones((n, 1)) / np.sqrt(n)
    return np.c_[np.kron(helmert(R), one(C)), np.kron(one(R), helmert(C)), np.kron(helmert(R), helmert(C))]


def ipf(x, w, m, it=5000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def neff(s):
    s = np.asarray(s, float) / np.sum(s); return 1.0 / np.sum(s ** 2)


def shares(n, balance, family):
    """n shares, largest first, with effective number = balance * n"""
    if balance >= 1 - 1e-12 or n == 1: return np.full(n, 1.0 / n)
    target = balance * n
    def make(a):                                        # a in (0, 1): 0 = equal, 1 = extreme
        if family == "geometric": s = (1 - a) ** np.arange(n, dtype=float)
        elif family == "dominant": s = np.r_[1.0, np.full(n - 1, (1 - a) / (1 + a * (n - 1)))]      # one large
        elif family == "tiny": s = np.r_[np.ones(n - 1), 1 - a]                                      # one small
        else: raise ValueError(family)
        return s / s.sum()
    lo, hi = 0.0, 1 - 1e-9
    if neff(make(hi)) > target: raise ValueError(f"family '{family}' cannot reach balance {balance} with {n} categories (minimum {neff(make(hi)) / n:.3f})")
    for _ in range(200):
        mid = (lo + hi) / 2
        if neff(make(mid)) > target: lo = mid
        else: hi = mid
    return np.sort(make((lo + hi) / 2))[::-1]


def make_design(rng=None, R=5, C=5, volume=("total", 3e4), eps=1.0, balance_row=1.0, balance_col=1.0, family_row="geometric",
                family_col="geometric", align="same", jitter=0.0, kappa=0.0, sigma_b=0.5, sigma_main=None, conflict=0.0):
    rng = rng or np.random.default_rng(0); D = R * C; nm = R + C - 2
    unit, val = volume
    N = {"total": val, "margin": val * R, "cell": val * R * C}[unit]
    w = shares(R, balance_row, family_row); m = shares(C, balance_col, family_col)
    if align == "reversed": m = m[::-1]
    elif align == "shifted": m = np.roll(m, 1)
    elif align != "same": raise ValueError(align)
    if jitter > 0: w = w * np.exp(rng.normal(0, jitter, R)); m = m * np.exp(rng.normal(0, jitter, C))
    w, m = N * w / w.sum(), N * m / m.sum()
    V = split_basis(R, C)
    centre = ipf(np.exp(kappa * np.eye(R, C)), w, m)
    s_main = sigma_b if sigma_main is None else sigma_main
    sig = np.r_[np.full(nm, s_main), np.full(D - 1 - nm, sigma_b)]
    if conflict > 0:
        u = rng.standard_normal(nm); u /= np.linalg.norm(u)
        centre = np.exp(np.log(centre) + (V[:, :nm] @ (conflict * s_main * u)).reshape(R, C)); centre *= N / centre.sum()
    lo = np.maximum(0, w[:, None] + m[None, :] - N); hi = np.minimum(w[:, None], m[None, :])
    feat = dict(R=R, C=C, total=N, eps=eps, margin_min=float(min(w.min(), m.min())), margin_max=float(max(w.max(), m.max())),
                row_margin_mean=N / R, col_margin_mean=N / C, cell_mean=N / D, cell_min_centre=float(centre.min()),
                neff_row=neff(w), neff_col=neff(m), balance_row=neff(w) / R, balance_col=neff(m) / C,
                interior_share=(R - 1) * (C - 1) / D, diag_mass=float(np.trace(centre[:min(R, C), :min(R, C)]) / N),
                frechet_width=float(np.mean((hi - lo) / hi)),
                tight_index_max=float(sigma_b ** 2 * np.sqrt(max(w.max(), m.max())) / eps),      # direct schemes: tightest margin
                tight_index_min=float(sigma_b ** 2 * np.sqrt(min(w.min(), m.min())) / eps),
                sigma_b=sigma_b, sigma_main=s_main, kappa=kappa, conflict=conflict)
    return dict(w=w, m=m, centre=centre, sigma=sig, V=V, N=N, eps=eps, features=feat)


if __name__ == "__main__":
    import pandas as pd
    pd.set_option("display.width", 250)
    show = ["total", "row_margin_mean", "cell_mean", "margin_min", "margin_max", "balance_row", "balance_col", "diag_mass", "frechet_width", "tight_index_max", "tight_index_min"]
    rows = []
    for lab, kw in [("baseline 5x5, total 30,000", {}), ("8x8, same total", dict(R=8, C=8)), ("8x8, same count per margin", dict(R=8, C=8, volume=("margin", 6000))),
                    ("8x8, same count per cell", dict(R=8, C=8, volume=("cell", 1200))), ("3x8, same count per cell", dict(R=3, C=8, volume=("cell", 1200))),
                    ("rows unbalanced (geometric)", dict(balance_row=0.5)), ("rows unbalanced (one dominant)", dict(balance_row=0.5, family_row="dominant")),
                    ("rows unbalanced (one tiny)", dict(balance_row=0.85, family_row="tiny")), ("both unbalanced, aligned", dict(balance_row=0.5, balance_col=0.5)),
                    ("both unbalanced, reversed", dict(balance_row=0.5, balance_col=0.5, align="reversed")),
                    ("diagonal centre", dict(kappa=3.0)), ("diagonal centre, both unbalanced, reversed", dict(kappa=3.0, balance_row=0.5, balance_col=0.5, align="reversed"))]:
        f = make_design(**kw)["features"]; rows.append(dict(design=lab, **{k: f[k] for k in show}))
    print(pd.DataFrame(rows).set_index("design").round(2).to_string())

"""Section 2 of the write-up: each allocation function is a one-to-one map from R^K onto the tables
with the given margins, with the stated Jacobian.

For every scheme and table size this checks, at random margins and parameters:
  margins   largest relative error in a row or column total
  positive  smallest cell (as a share of the table total)
  inverse   largest |inv(alloc(lam)) - lam|
  jacobian  largest |stated log Jacobian - finite-difference log Jacobian|
  onto      largest relative error in alloc(inv(T)) - T over random interior tables T
            (and the share of tables that cannot be reached, for the smoothed bounds)

Run from the repository root:  python3 experiments/01_allocation_works.py
Needs cmdstanpy and CmdStan. Writes experiments/output/01_allocation_works.md
"""
import os, sys, numpy as np, logging
import cmdstanpy

logging.getLogger("cmdstanpy").disabled = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
model = cmdstanpy.CmdStanModel(stan_file=os.path.join(ROOT, "stan", "check_alloc.stan"),
                               stanc_options={"include-paths": [os.path.join(ROOT, "stan")]})

SCHEMES = [("between bounds, sharp", 1, 0.0), ("between bounds, smoothed (delta 0.05)", 1, 0.05),
           ("cell logit", 2, 0.0), ("adjusted row", 3, 0.0), ("adjusted table", 4, 0.0)]
SIZES = [(2, 2), (3, 3), (3, 5), (5, 3), (5, 5), (7, 7)]
P, P2, H, TOTAL = 12, 200, 1e-4, 1000.0


def run(R, C, scheme, delta, row_order, rem_col, w, m, lam, Tin):
    data = dict(R=R, C=C, scheme=scheme, row_order=row_order, rem_col=rem_col, delta=delta,
                N=len(lam), w=w, m=m, lam=lam, N2=len(Tin), Tin=Tin)
    fit = model.sample(data=data, fixed_param=True, iter_sampling=1, chains=1, sig_figs=18, show_progress=False)
    g = lambda name: fit.stan_variable(name)[0]
    return g("T"), g("lj"), g("lam_back"), g("lam_in"), g("T_back")


def check(R, C, scheme, delta, row_order, rem_col, rng):
    K = (R - 1) * (C - 1)
    w0 = rng.dirichlet(np.full(R, 2.0), P) * TOTAL
    m0 = rng.dirichlet(np.full(C, 2.0), P) * TOTAL
    lam0 = rng.normal(0, 1.5, (P, K))
    # each base point followed by its 2K finite-difference neighbours
    w, m, lam = [], [], []
    for p in range(P):
        w.append(w0[p]); m.append(m0[p]); lam.append(lam0[p])
        for k in range(K):
            for s in (1, -1):
                l = lam0[p].copy(); l[k] += s * H
                w.append(w0[p]); m.append(m0[p]); lam.append(l)
    Tin = rng.dirichlet(np.ones(R * C), P2).reshape(P2, R, C) * TOTAL
    T, lj, lam_back, lam_in, T_back = run(R, C, scheme, delta, row_order, rem_col, w, m, lam, Tin)
    step = 2 * K + 1
    out = dict(margins=0.0, positive=np.inf, inverse=0.0, jacobian=0.0)
    for p in range(P):
        b = p * step
        Tb = T[b]
        out["margins"] = max(out["margins"], np.abs(Tb.sum(1) - w0[p]).max() / TOTAL, np.abs(Tb.sum(0) - m0[p]).max() / TOTAL)
        out["positive"] = min(out["positive"], Tb.min() / TOTAL)
        out["inverse"] = max(out["inverse"], np.abs(lam_back[b] - lam0[p]).max())
        J = np.zeros((K, K))
        for k in range(K):
            J[:, k] = (T[b + 1 + 2 * k][:R - 1, :C - 1].ravel() - T[b + 2 + 2 * k][:R - 1, :C - 1].ravel()) / (2 * H)
        out["jacobian"] = max(out["jacobian"], abs(np.linalg.slogdet(J)[1] - lj[b]))
    ok = np.isfinite(lam_in).all(1)
    out["unreachable"] = 1 - ok.mean()
    out["onto"] = (np.abs(T_back[ok] - Tin[ok]).max() / TOTAL) if ok.any() else np.nan
    return out


def main():
    rng = np.random.default_rng(20261005)
    lines = ["# Check that each allocation function works", "",
             f"{P} random margins and parameter vectors per row (parameters ~ N(0, 1.5)), table total {TOTAL:g}; "
             f"{P2} random interior tables for the onto check. Orders: natural = rows in order, last row and last "
             "column by subtraction; shuffled = a random row order and a random subtraction column in each row.", "",
             "| Scheme | Size | Order | Margins | Smallest cell | Inverse | Jacobian | Onto | Unreachable |",
             "|---|---|---|---|---|---|---|---|---|"]
    worst = {}
    for name, scheme, delta in SCHEMES:
        for R, C in SIZES:
            orders = [("natural", list(range(1, R + 1)), [C] * R)]
            if scheme != 4 and (R, C) != (2, 2):
                orders.append(("shuffled", [int(x) for x in rng.permutation(R) + 1], [int(x) for x in rng.integers(1, C + 1, R)]))
            for oname, ro, rc in orders:
                o = check(R, C, scheme, delta, ro, rc, rng)
                lines.append(f"| {name} | {R}x{C} | {oname} | {o['margins']:.1e} | {o['positive']:.1e} | {o['inverse']:.1e} | "
                             f"{o['jacobian']:.1e} | {o['onto']:.1e} | {o['unreachable']:.1%} |")
                print(lines[-1], flush=True)
                for k in ("margins", "inverse", "jacobian", "onto", "unreachable"):
                    worst[(name, k)] = max(worst.get((name, k), 0.0), o[k])
                worst[(name, "positive")] = min(worst.get((name, "positive"), np.inf), o["positive"])
    lines += ["", "## Worst case by scheme", "", "| Scheme | Margins | Smallest cell | Inverse | Jacobian | Onto | Unreachable |", "|---|---|---|---|---|---|---|"]
    for name, _, _ in SCHEMES:
        lines.append(f"| {name} | {worst[(name, 'margins')]:.1e} | {worst[(name, 'positive')]:.1e} | {worst[(name, 'inverse')]:.1e} | "
                     f"{worst[(name, 'jacobian')]:.1e} | {worst[(name, 'onto')]:.1e} | {worst[(name, 'unreachable')]:.1%} |")
        print(lines[-1])
    with open(os.path.join(ROOT, "experiments", "output", "01_allocation_works.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

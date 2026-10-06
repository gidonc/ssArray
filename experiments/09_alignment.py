"""Is the trouble with the direct parameterisations linear (a matter of alignment) or nonlinear (curvature)?

With an isotropic prior the density on tables is the same for every orthonormal ILR basis V, so V only changes how the
coordinates line up with the stiff margin directions.  A dense metric removes any linear misalignment; a basis whose
leading coordinates span the margin directions (at the prior centre) removes it by construction.  Whatever difficulty is
left under either is curvature.
Variants on the Scottish grand total, prior centred on independence:
  ILR Helmert / ILR aligned / log cells   x   diag / dense metric,   with position logit (scaled margins) as reference.
Run:  python3 experiments/09_alignment.py     Writes experiments/output/09_alignment.csv and .md
"""
import os, time, importlib.util, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2


def aligned_basis(t):
    """orthonormal basis of the clr subspace; the first R+C-2 columns span the margin directions d log margins / d log cells at t"""
    R, C = t.shape; D = R * C; a = []
    for r in range(R):
        v = np.zeros((R, C)); v[r] = t[r] / t[r].sum(); a.append(v.ravel())
    for c in range(C):
        v = np.zeros((R, C)); v[:, c] = t[:, c] / t[:, c].sum(); a.append(v.ravel())
    a = np.array(a); a = a - a.mean(1, keepdims=True)
    U, s, _ = np.linalg.svd(np.c_[a.T, np.eye(D) - 1.0 / D], full_matrices=False)     # margin directions first, then the rest of clr
    k = R + C - 2
    Q, _ = np.linalg.qr(np.c_[np.linalg.svd(a.T, full_matrices=False)[0][:, :k], (np.eye(D) - 1.0 / D)])
    Q = Q[:, :D]
    keep = [j for j in range(Q.shape[1]) if abs(Q[:, j].sum()) < 1e-8][:D - 1]
    return Q[:, keep]


def main():
    T = r4.load(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R, C = T.shape; D = R * C
    T0 = np.outer(T.sum(1), T.sum(0)) / T.sum()
    Va = aligned_basis(T0)
    assert Va.shape == (D, D - 1) and np.allclose(Va.T @ Va, np.eye(D - 1), atol=1e-8) and np.allclose(Va.sum(0), 0, atol=1e-8)
    variants = [("ILR Helmert", 0, None), ("ILR aligned", 0, Va), ("log cells", 5, None), ("position logit", 1, None)]
    rows = []; out = os.path.join(ROOT, "experiments", "output", "09_alignment.csv")
    for eps in (1.0, 0.1):
        for sb in (1.0, 2.0):
            for name, p, Vx in variants:
                for metric in ("diag_e", "dense_e"):
                    if p == 1 and metric == "dense_e": continue
                    for seed in (1, 2):
                        d = e2.base_data(T, eps, 1 if p == 1 else 0)
                        if Vx is not None: d["V"] = Vx
                        d.update(param=p, mu_b=d["V"].T @ np.log(T0).ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                        row = dict(eps=eps, sigma_b=sb, variant=name, metric=metric, seed=seed); t0 = time.time()
                        try:
                            f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=seed, metric=metric,
                                                inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False, timeout=200)
                            s = f.summary(); c = s[s.index.str.startswith("log_T[")]; dr = f.draws_pd()
                            row.update(seconds=time.time() - t0, max_rhat=float(c["R_hat"].max()), min_ess=float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()),
                                       div=int(dr["divergent__"].sum()), leapfrogs=float(dr["n_leapfrog__"].mean()))
                        except Exception as ex:
                            row.update(seconds=time.time() - t0, max_rhat=np.inf, min_ess=0.0, error=type(ex).__name__)
                        rows.append(row); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
                        pd.DataFrame(rows).to_csv(out, index=False)
    df = pd.DataFrame(rows); df["fail"] = df.max_rhat > 1.05; df["ess_per_s"] = df.min_ess / df.seconds
    g = df.groupby(["eps", "sigma_b", "variant", "metric"]).agg(failed=("fail", "sum"), runs=("fail", "size"), max_rhat=("max_rhat", "max"),
            min_ess=("min_ess", "median"), leapfrogs=("leapfrogs", "median"), seconds=("seconds", "median"), ess_per_s=("ess_per_s", "median")).round(2)
    open(out.replace(".csv", ".md"), "w").write("# Alignment test, Scottish grand total\n\n" + g.to_markdown() + "\n")
    print(g.to_string())


if __name__ == "__main__":
    main()

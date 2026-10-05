"""Curvature of the six coordinate systems on the Scottish grand total.

All six share one posterior density, so any difference is the geometry of the coordinates. For each, take posterior
draws theta_i, compute the Hessian H_i of the log density in theta (BridgeStan), and look at it after whitening by the
posterior covariance M = cov(theta) (the best a single dense metric can do): G_i = -L' H_i L, M = L L'.
  cond     median over draws of the condition number of G_i (|eigenvalues|)
  vary     median over draws of ||G_i - mean(G)||_F / ||mean(G)||_F : how much the curvature moves around the posterior,
           which a constant metric cannot absorb (0 = a Gaussian in these coordinates)
  neg      share of draws where G_i has a negative eigenvalue (non-log-concave there)
Run from the repository root:  python3 experiments/06_curvature.py [--eps 1] [--sigma 1,2] [--params 0,5,1,2,3,4] [--draws 40]
Needs BridgeStan (pip install bridgestan, plus a clone of github.com/roualdes/bridgestan with --recurse-submodules,
path in BRIDGESTAN).  Writes experiments/output/06_curvature.csv and .md.
"""
import os, sys, json, argparse, importlib.util, numpy as np, pandas as pd
import bridgestan as bs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4)
e2 = r4.e2

bs.set_bridgestan_path(os.environ["BRIDGESTAN"])
SO = os.path.join(ROOT, "stan", "single_area_model.so")
if not os.path.exists(SO):
    bs.compile_model(os.path.join(ROOT, "stan", "single_area.stan"), stanc_args=[f"--include-paths={ROOT}/stan"])


def jsonable(d):
    return json.dumps({k: (np.asarray(v).tolist() if not np.isscalar(v) else (float(v) if isinstance(v, (float, np.floating)) else int(v))) for k, v in d.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", default="data/scot_2007_all_areas.csv")
    ap.add_argument("--eps", default="1"); ap.add_argument("--sigma", default="1,2")
    ap.add_argument("--params", default="0,5,1,2,3,4"); ap.add_argument("--draws", type=int, default=40)
    ap.add_argument("--label", default="curvature")
    a = ap.parse_args()
    T = r4.load(os.path.join(ROOT, a.path)); R, C = T.shape; D = R * C
    rows = []
    for eps in [float(x) for x in a.eps.split(",")]:
        for sb in [float(x) for x in a.sigma.split(",")]:
            for p in [int(x) for x in a.params.split(",")]:
                sc = 0 if p in (0, 5) else 1
                d = e2.base_data(T, eps, sc)
                d.update(param=p, mu_b=d["V"].T @ r4.centre_table(T, "indep").ravel(), sigma_b=np.full(D - 1, sb), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=1,
                                    inits=[{"theta": e2.start(p, d)}] * 2, show_progress=False, show_console=False)
                s = f.summary(); c = s[s.index.str.startswith("log_T[")]
                rhat = float(c["R_hat"].max())
                th = f.stan_variable("theta")
                if th.ndim == 3: th = th.reshape(-1, th.shape[-1])
                model = bs.StanModel(SO, jsonable({k: v for k, v in d.items()}), seed=1)
                M = np.cov(th.T); L = np.linalg.cholesky(M + 1e-12 * np.trace(M) / D * np.eye(D))
                idx = np.linspace(0, len(th) - 1, a.draws).astype(int)
                Gs = []
                for i in idx:
                    _, _, H = model.log_density_hessian(np.ascontiguousarray(th[i], dtype=np.float64), jacobian=True)
                    Gs.append(-(L.T @ H @ L))
                Gs = np.array(Gs); Gs = (Gs + Gs.transpose(0, 2, 1)) / 2
                ev = np.linalg.eigvalsh(Gs)
                cond = np.median(np.abs(ev).max(1) / np.abs(ev).min(1))
                Gm = Gs.mean(0)
                vary = np.median([np.linalg.norm(G - Gm) / np.linalg.norm(Gm) for G in Gs])
                neg = float((ev.min(1) < 0).mean())
                rows.append(dict(eps=eps, sigma_b=sb, param=r4.NAMES[p], rhat=rhat, min_ess=float(c["ESS_bulk" if "ESS_bulk" in c else "N_Eff"].min()),
                                 cond=cond, vary=vary, neg=neg, max_eig=float(np.median(ev.max(1))), min_eig=float(np.median(ev.min(1)))))
                print(rows[-1], flush=True)
                pd.DataFrame(rows).to_csv(os.path.join(ROOT, "experiments", "output", f"06_{a.label}.csv"), index=False)
    df = pd.DataFrame(rows)
    with open(os.path.join(ROOT, "experiments", "output", f"06_{a.label}.md"), "w") as fh:
        fh.write("# Curvature after whitening by the posterior covariance\n\nvary = how much the Hessian moves around the posterior (0 = Gaussian); "
                 "cond = condition number of the whitened Hessian; neg = share of draws with a negative eigenvalue.\n"
                 "Only rows with rhat < 1.05 describe a converged posterior.\n\n" + df.round(3).to_markdown(index=False) + "\n")
    print(df.round(3).to_string())


if __name__ == "__main__":
    main()

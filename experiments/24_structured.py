"""Large tables after three changes to the model:
  split_basis = 1       the prior on the split basis evaluated without the dense RC x (RC-1) matrix
  centred_interior = 1  adjusted row and adjusted table in contrasts with no reference cell
  adjusted table        Jacobian from an (R+C-1) square determinant instead of an (R-1)(C-1) square one
Part 'grad': time per gradient by size, dense basis against split basis, all six parameterisations (short runs).
Part 'big' : as experiment 23 (2 chains, 500 + 500 draws, limit 900 s) with the split basis and centred interior, adding 48x48.
Designs as in 23: about 6,000 per row margin, margins jittered, independence centre, sigma_b 0.5, eps 1, largest-column
margin coordinates.
Run:  python3 experiments/24_structured.py grad|big     Writes experiments/output/24_structured_<part>.csv and .md
"""
import os, sys, time, importlib.util, numpy as np, pandas as pd, arviz as az
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "experiments"))
import simdesign as sd, margincoords as mc
spec = importlib.util.spec_from_file_location("r4", os.path.join(ROOT, "experiments", "04_real_table.py"))
r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4); e2 = r4.e2


def fit(R, C, p, split, centred, warm, draws, limit):
    d = sd.make_design(np.random.default_rng(R * 100 + C), R=R, C=C, volume=("margin", 6000), jitter=0.15); w, m = d["w"], d["m"]; D = R * C
    T = np.outer(w, m) / d["N"]; K, sc = mc.margin_K(w, m, 1.0, "largest", 1)
    dd = e2.base_data(T, 1.0, 0 if p in (0, 5) else sc); dd["K_margin"] = K; dd["split_basis"] = split; dd["centred_interior"] = centred
    dd["V"] = np.zeros((0, 0)) if split else d["V"]
    dd.update(param=p, mu_b=d["V"].T @ np.log(T).ravel(), sigma_b=d["sigma"], mu_logv=float(np.log(d["N"])), sigma_logv=1.0)
    row = dict(R=R, C=C, cells=D, param=r4.NAMES[p], basis="split" if split else "dense", centred_interior=centred); t0 = time.time()
    try:
        f = e2.model.sample(data=dd, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=draws, seed=1, save_warmup=True,
                            inits=[{"theta": e2.start(p, dd)}] * 2, show_progress=False, show_console=False, timeout=limit)
        wall = time.time() - t0
        dr = f.draws_pd(inc_warmup=True); post = dr[dr["iter__"] > warm]
        cols = [c for c in dr.columns if c.startswith("log_T[")]; step = max(1, len(cols) // 200)
        x = np.stack([post[post["chain__"] == ch][cols].values for ch in sorted(post["chain__"].unique())])
        ess = float(np.min([az.ess(x[:, :, j], method="bulk") for j in range(0, x.shape[2], step)]))
        rh = float(np.max([az.rhat(x[:, :, j]) for j in range(0, x.shape[2], step)]))
        grads = dr.groupby("chain__")["n_leapfrog__"].sum().mean()
        row.update(seconds=wall, leapfrogs=float(post["n_leapfrog__"].mean()), ms_per_gradient=1000 * wall / grads, min_ess=ess, max_rhat=rh, ess_per_s=ess / wall,
                   div=int(post["divergent__"].sum()), status="ok" if rh < 1.05 else "not converged")
    except Exception as ex:
        row.update(seconds=time.time() - t0, status=type(ex).__name__)
    return row


def main():
    part = sys.argv[1]; rows = []; out = os.path.join(ROOT, "experiments", "output", "24_structured_" + part)
    if part == "grad":
        for R in (3, 5, 7, 10, 16, 24, 32):
            for p in (2, 1, 3, 4, 0, 5):
                for split in (0, 1):
                    rows.append(fit(R, R, p, split, 1, 150, 150, 600)); pd.DataFrame(rows).to_csv(out + ".csv", index=False)
                    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        df = pd.DataFrame(rows); df["size"] = df.R.astype(str) + "x" + df.C.astype(str)
        g = df.pivot_table(index=["R", "size"], columns=["param", "basis"], values="ms_per_gradient").round(3)
        open(out + ".md", "w").write("# Time per gradient (ms): dense basis against split basis\n\n" + g.to_markdown() + "\n")
    else:
        for R, C in [(10, 10), (5, 40), (16, 16), (24, 24), (32, 32), (48, 48)]:
            for p in (2, 1, 3, 4, 5, 0):
                rows.append(fit(R, C, p, 1, 1, 500, 500, 900)); pd.DataFrame(rows).to_csv(out + ".csv", index=False)
                print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        df = pd.DataFrame(rows)
        open(out + ".md", "w").write("# Large tables with the split basis and centred interior: 2 chains, 500 + 500 draws, limit 900 s\n\n" + df.round(2).to_markdown(index=False) + "\n")


if __name__ == "__main__":
    main()

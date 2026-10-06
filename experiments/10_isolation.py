"""Are the tight directions isolated in a few parameters, or spread over many?  One table, one density, six coordinate systems.

At a table t the density has a precision in log-cell space, P_x = prior (ILR, log total) + margin penalty.  In a
parameterisation theta it becomes P = J' P_x J, J = d log cells / d theta (curvature terms of the map left out: this is the
linear picture at t).  The k = R + C - 1 largest eigenvalues of P are the directions the margins make tight.  For each:
  cond        condition number of P
  cond_diag   the same after rescaling each parameter (unit diagonal): what a diagonal metric is left with
  carriers    median number of parameters carrying a tight direction (participation ratio; 1 = one parameter each)
  in_k        share of the tight subspace lying in the k most involved parameters (1 = isolated in k parameters)
  rotation    mean principal angle (degrees) between the tight subspace at t and at tables one prior sd away with the same
              margins: 0 = the tight directions stay in the same parameters across the posterior
Tables: the independence table and the actual table of the Scottish grand total (same margins).
Run:  BRIDGESTAN=/path/to/bridgestan python3 experiments/10_isolation.py     Writes experiments/output/10_isolation.csv and .md
"""
import os, importlib.util, numpy as np, pandas as pd
from scipy.optimize import least_squares
import bridgestan as bs
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("c6", os.path.join(ROOT, "experiments", "06_curvature.py"))
c6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c6)
r4, e2, jsonable, SO = c6.r4, c6.e2, c6.jsonable, c6.SO
EPS, SB = 0.1, 1.0


def ipf(x, w, m, it=5000):
    x = x.copy()
    for _ in range(it):
        x *= (w / x.sum(1))[:, None]; x *= (m / x.sum(0))[None, :]
    return x


def P_x(t, V, obs_w, obs_m):
    R, C = t.shape; D = R * C; a = []
    for r in range(R):
        v = np.zeros((R, C)); v[r] = t[r]; a.append(v.ravel())                    # d row margin / d log cells
    for c in range(C):
        v = np.zeros((R, C)); v[:, c] = t[:, c]; a.append(v.ravel())
    a = np.array(a); W = np.diag(1 / (EPS ** 2 * np.r_[obs_w, obs_m]))
    s = (t / t.sum()).ravel()[:, None]                                            # d log total / d log cells
    return V @ V.T / SB ** 2 + s @ s.T + a.T @ W @ a


def main():
    T = r4.load(os.path.join(ROOT, "data", "scot_2007_all_areas.csv")); R, C = T.shape; D = R * C; k = R + C - 1
    w, m = T.sum(1), T.sum(0); Ti = np.outer(w, m) / T.sum()
    rng = np.random.default_rng(1); rows = []
    for cname, t0 in (("independence", Ti), ("actual", T)):
        V = e2.helmert(D)
        near = [ipf(np.exp(np.log(t0).ravel() + V @ rng.normal(0, SB, D - 1)).reshape(R, C), w, m) for _ in range(4)]
        for p in (0, 5, 1, 2, 3, 4):
            for sc in ((0,) if p in (0, 5) else (1, 0)):
                d = e2.base_data(T, EPS, sc)
                d.update(param=p, mu_b=V.T @ np.log(t0).ravel(), sigma_b=np.full(D - 1, SB), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
                model = bs.StanModel(SO, jsonable(d), seed=1)
                logT = lambda th: np.log(model.param_constrain(np.ascontiguousarray(th, dtype=np.float64), include_tp=True)[D:2 * D].reshape((R, C), order="F")).ravel()

                def theta_of(t):
                    lt = np.log(t).ravel()
                    if p == 0: return np.r_[V.T @ lt, np.log(t.sum())]
                    if p == 5: return lt
                    f = lambda lam: logT(np.r_[np.zeros(k), lam]) - lt             # margins are the observed ones: margin parameters 0
                    s = least_squares(f, np.zeros(D - k), xtol=1e-14, ftol=1e-14, gtol=1e-14)
                    assert np.abs(s.fun).max() < 1e-6, (r4.NAMES[p], np.abs(s.fun).max())
                    return np.r_[np.zeros(k), s.x]

                def tight(t):
                    th = theta_of(t); h = 1e-5
                    J = np.array([(logT(th + h * e) - logT(th - h * e)) / (2 * h) for e in np.eye(D)]).T
                    P = J.T @ P_x(t, V, w, m) @ J; P = (P + P.T) / 2
                    ev, U = np.linalg.eigh(P)
                    return P, ev, U[:, -k:]

                P, ev, U = tight(t0)
                dg = np.sqrt(np.diag(P)); evd = np.linalg.eigvalsh(P / np.outer(dg, dg))
                pr = [(u ** 2).sum() ** 2 / (u ** 4).sum() for u in U.T]
                load = (U ** 2).sum(1)                                              # each parameter's share of the tight subspace (sums to k)
                rot = [np.degrees(np.arccos(np.clip(np.linalg.svd(U.T @ tight(t)[2], compute_uv=False), -1, 1))).mean() for t in near]
                rows.append(dict(table=cname, param=r4.NAMES[p], scale_margins=sc, cond=ev[-1] / ev[0], cond_diag=evd[-1] / evd[0],
                                 carriers=float(np.median(pr)), carriers_max=float(np.max(pr)), in_k=float(np.sort(load)[-k:].sum() / k),
                                 rotation=float(np.mean(rot))))
                print(rows[-1], flush=True); del model
    df = pd.DataFrame(rows); out = os.path.join(ROOT, "experiments", "output", "10_isolation")
    df.to_csv(out + ".csv", index=False)
    f = df.copy(); f["cond"] = f.cond.map("{:.1e}".format); f["cond_diag"] = f.cond_diag.map("{:.1e}".format)
    open(out + ".md", "w").write(f"# Where the tight directions sit, Scottish grand total (eps {EPS}, sigma_b {SB}, k = {k} tight directions of {D})\n\n" + f.round(2).to_markdown(index=False) + "\n")
    print(f.round(2).to_string())


if __name__ == "__main__":
    main()

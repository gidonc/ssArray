"""The single-area model in six coordinate systems: check the maps, then check they sample the same posterior.

Part 1 (maps): for each parameterisation, the stated log Jacobian of (parameters -> cells) against finite
differences, on the 3x3 and 7x7 versions of the example table.

Part 2 (sampling): the same density, sampled in each parameterisation at a loose and a tight margin
penalty. Reports sampler behaviour and how far each posterior is from the adjusted-row one.
Short runs: this is a check that the models work, not the curvature comparison.

Example table: one Scottish Parliament constituency (data/scot_area60.csv), party-list vote (rows)
by constituency vote (columns); the 3x3 version collapses to SNP, Labour, others.

Run from the repository root:  python3 experiments/02_single_area_models.py
Writes experiments/output/02_single_area_models.md
"""
import os, csv, time, logging, numpy as np
import cmdstanpy

logging.getLogger("cmdstanpy").disabled = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAN = os.path.join(ROOT, "stan")
opts = {"include-paths": [STAN]}
check = cmdstanpy.CmdStanModel(stan_file=os.path.join(STAN, "check_single_area.stan"), stanc_options=opts)
model = cmdstanpy.CmdStanModel(stan_file=os.path.join(STAN, "single_area.stan"), stanc_options=opts)

PARAMS = [(0, "ILR + log volume"), (5, "log cells"), (1, "between bounds"), (2, "cell logit"),
          (3, "adjusted row"), (4, "adjusted table")]
REF = 3


def load_table():
    T = np.zeros((7, 7))
    with open(os.path.join(ROOT, "data", "scot_area60.csv")) as f:
        for row in csv.DictReader(f):
            T[int(row["row_no"]) - 1, int(row["col_no"]) - 1] = float(row["votes"])
    return T


def collapse(T, groups):
    return np.array([[T[np.ix_(a, b)].sum() for b in groups] for a in groups])


def helmert(D):
    V = np.zeros((D, D - 1))
    for k in range(1, D):
        V[:k, k - 1] = 1 / np.sqrt(k * (k + 1.0))
        V[k, k - 1] = -k / np.sqrt(k * (k + 1.0))
    return V


def base_data(T, eps):
    R, C = T.shape
    w, m = T.sum(1), T.sum(0)
    order = [int(i) + 1 for i in np.argsort(w)]                    # smallest row first, largest by subtraction
    return dict(R=R, C=C, obs_w=w, obs_m=m, eps=eps, V=helmert(R * C), row_order=order,
                rem_col=[r + 1 if r < C else C for r in range(R)], delta=0.0)


def start(param, d):
    """a starting point on the margins: the independence table"""
    R, C = d["R"], d["C"]
    T0 = np.outer(d["obs_w"], d["obs_m"]) / d["obs_w"].sum()
    if param == 0:
        return np.r_[d["V"].T @ np.log(T0).ravel(), np.log(T0.sum())]
    if param == 5:
        return np.log(T0).ravel()
    return np.zeros(R * C)


def jacobian_check(T, eps, rng, P=6, H=1e-5):
    out = []
    d = base_data(T, eps)
    D = d["R"] * d["C"]
    for param, name in PARAMS:
        th0 = [start(param, d) + rng.normal(0, 0.5, D) for _ in range(P)]
        thetas = []
        for t in th0:
            thetas.append(t)
            for k in range(D):
                for s in (1, -1):
                    u = t.copy(); u[k] += s * H; thetas.append(u)
        fit = check.sample(data=dict(d, param=param, N=len(thetas), theta=thetas), fixed_param=True,
                           iter_sampling=1, chains=1, sig_figs=18, show_progress=False)
        Tt, lj = fit.stan_variable("T")[0], fit.stan_variable("lj")[0]
        diff = []
        for p in range(P):
            b = p * (2 * D + 1)
            Jm = np.stack([(Tt[b + 1 + 2 * k].ravel() - Tt[b + 2 + 2 * k].ravel()) / (2 * H) for k in range(D)], 1)
            diff.append(np.linalg.slogdet(Jm)[1] - lj[b])
        diff = np.array(diff)
        out.append((name, diff.mean(), np.abs(diff - diff.mean()).max()))
    return out


def sample(T, eps, param, seed=1, warm=500, draws=500):
    d = base_data(T, eps)
    D = d["R"] * d["C"]
    d.update(param=param, mu_b=np.zeros(D - 1), sigma_b=np.full(D - 1, 2.0), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
    t0 = time.time()
    fit = model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=draws, seed=seed,
                       inits=[{"theta": start(param, d)}] * 2, show_progress=False, show_console=False)
    secs = time.time() - t0
    s = fit.summary()
    cells = s[s.index.str.startswith("log_T[")]
    ess_col = "N_Eff" if "N_Eff" in s.columns else "ESS_bulk"
    dr = fit.draws_pd()
    leap = dr["n_leapfrog__"].mean()
    return dict(mean=fit.stan_variable("log_T").mean(0), sd=fit.stan_variable("log_T").std(0),
                step=dr["stepsize__"].mean(), leap=leap, div=int(dr["divergent__"].sum()),
                tree=float((dr["treedepth__"] >= 10).mean()), ess=float(cells[ess_col].min()),
                rhat=float(cells["R_hat"].max()), secs=secs, grads=leap * 2 * draws)


def main():
    rng = np.random.default_rng(20261005)
    T7 = load_table()
    T3 = collapse(T7, [[5], [1], [0, 2, 3, 4, 6]])                 # SNP, Labour, others
    lines = ["# The single-area model in six coordinate systems", "",
             "## Part 1: the stated Jacobian against finite differences", "",
             "Six random points per row. 'Offset' is the mean of (finite difference - stated); it should be zero, or a "
             "constant for the ILR parameterisation, whose Jacobian is stated up to a constant. 'Spread' is the largest "
             "deviation from that mean and should be near zero.", "",
             "| Size | Parameterisation | Offset | Spread |", "|---|---|---|---|"]
    for label, T in (("3x3", T3), ("7x7", T7)):
        for name, off, spread in jacobian_check(T, 0.3, rng):
            lines.append(f"| {label} | {name} | {off:.2e} | {spread:.1e} |")
            print(lines[-1], flush=True)
    lines += ["", "## Part 2: sampling the same posterior", "",
              "2 chains, 500 warm-up + 500 draws each, started at the independence table. Placeholder prior: ILR "
              "coordinates ~ normal(0, 2). Summaries are of the log cells. 'Gap' is the largest difference in a posterior mean from the "
              "adjusted-row run, in posterior sds; 'sd ratio' is the range of the ratio of posterior sds to that run. "
              "ESS is the smallest over cells, of 1,000 draws.", "",
              "| Size | eps | Parameterisation | Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | Seconds | Gap | sd ratio |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for label, T in (("3x3", T3), ("7x7", T7)):
        for eps in (1.0, 0.1):
            res = {p: sample(T, eps, p) for p, _ in PARAMS}
            ref = res[REF]
            for p, name in PARAMS:
                r = res[p]
                gap = np.abs(r["mean"] - ref["mean"]).__truediv__(ref["sd"]).max()
                ratio = r["sd"] / ref["sd"]
                lines.append(f"| {label} | {eps:g} | {name} | {r['step']:.3g} | {r['leap']:.0f} | {r['div']} | {r['tree']:.0%} | "
                             f"{r['ess']:.0f} | {r['rhat']:.3f} | {1000 * r['ess'] / r['grads']:.2f} | {r['secs']:.0f} | "
                             f"{gap:.2f} | {ratio.min():.2f}-{ratio.max():.2f} |")
                print(lines[-1], flush=True)
    lines += ["", "## Part 3: longer runs of the sequential parameterisations at 7x7", "",
              "2 chains, 1,000 warm-up + 3,000 draws each, eps = 1. Same columns; ESS is of 6,000 draws; "
              "gap and sd ratio are against the adjusted-table run.", "",
              "| Parameterisation | Step size | Leapfrogs | Divergences | Min ESS | Max Rhat | ESS per 1,000 gradients | Seconds | Gap | sd ratio |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    res = {p: sample(T7, 1.0, p, warm=1000, draws=3000) for p in (1, 2, 3, 4)}
    for p, name in PARAMS:
        if p not in res: continue
        r = res[p]
        gap = (np.abs(r["mean"] - res[4]["mean"]) / res[4]["sd"]).max()
        ratio = r["sd"] / res[4]["sd"]
        lines.append(f"| {name} | {r['step']:.3g} | {r['leap']:.0f} | {r['div']} | {r['ess']:.0f} | {r['rhat']:.3f} | "
                     f"{1000 * r['ess'] / r['grads']:.2f} | {r['secs']:.0f} | {gap:.2f} | {ratio.min():.2f}-{ratio.max():.2f} |")
        print(lines[-1], flush=True)
    with open(os.path.join(ROOT, "experiments", "output", "02_single_area_models.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

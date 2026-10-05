"""The single-area model in six coordinate systems: check the maps, check they sample the same posterior,
and record what each costs.

Part 1  the stated log Jacobian of (parameters -> cells) against finite differences.
Part 2  the same density sampled in each parameterisation, at a loose and a tight margin penalty,
        on 3x3, 5x5 and 7x7 versions of the example table. Efficiency per gradient and per second.
Part 3  longer 7x7 runs of the sequential parameterisations: do they agree?
Part 4  does scaling the margin parameters by the penalty width matter?
Part 5  cost of one gradient by table size.

Short runs on one table with a placeholder prior: a check that the models work, not the curvature comparison.

Example table: one Scottish Parliament constituency (data/scot_area60.csv), party-list vote (rows) by
constituency vote (columns). 3x3 = SNP, Labour, others; 5x5 = SNP, Labour, Conservative, Lib Dem, others.

Run from the repository root:  python3 experiments/02_single_area_models.py
Writes experiments/output/02_single_area_models.md
"""
import os, csv, logging, numpy as np
import cmdstanpy

logging.getLogger("cmdstanpy").disabled = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAN = os.path.join(ROOT, "stan")
opts = {"include-paths": [STAN]}
check = cmdstanpy.CmdStanModel(stan_file=os.path.join(STAN, "check_single_area.stan"), stanc_options=opts)
model = cmdstanpy.CmdStanModel(stan_file=os.path.join(STAN, "single_area.stan"), stanc_options=opts)

PARAMS = [(0, "ILR + log volume"), (5, "log cells"), (1, "between bounds"), (2, "cell logit"),
          (3, "adjusted row"), (4, "adjusted table")]
NAME = dict(PARAMS)


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


def base_data(T, eps, scale):
    R, C = T.shape
    w, m = T.sum(1), T.sum(0)
    order = [int(i) + 1 for i in np.argsort(w)]                    # smallest row first, largest by subtraction
    return dict(R=R, C=C, obs_w=w, obs_m=m, eps=eps, scale_margins=scale, V=helmert(R * C), row_order=order,
                rem_col=[r + 1 if r < C else C for r in range(R)], delta=0.0)


def start(param, d):
    """a starting point on the margins: the independence table"""
    T0 = np.outer(d["obs_w"], d["obs_m"]) / d["obs_w"].sum()
    if param == 0:
        return np.r_[d["V"].T @ np.log(T0).ravel(), np.log(T0.sum())]
    if param == 5:
        return np.log(T0).ravel()
    return np.zeros(d["R"] * d["C"])


def jacobian_check(T, eps, scale, rng, P=6, H=1e-5):
    out = []
    d = base_data(T, eps, scale)
    D = d["R"] * d["C"]
    for param, name in PARAMS:
        if scale == 1 and param in (0, 5):
            continue                                               # the switch does nothing for these
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


def sample(T, eps, param, scale=0, seed=1, warm=500, draws=500):
    d = base_data(T, eps, scale)
    D = d["R"] * d["C"]
    d.update(param=param, mu_b=np.zeros(D - 1), sigma_b=np.full(D - 1, 2.0), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
    fit = model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=warm, iter_sampling=draws, seed=seed,
                       inits=[{"theta": start(param, d)}] * 2, show_progress=False, show_console=False)
    s = fit.summary()
    cells = s[s.index.str.startswith("log_T[")]
    ess = float(cells["N_Eff" if "N_Eff" in s.columns else "ESS_bulk"].min())
    dr = fit.draws_pd()
    grads = float(dr["n_leapfrog__"].sum())                        # gradients in the sampling phase, both chains
    t_warm = sum(t["warmup"] for t in fit.time)
    t_samp = sum(t["sampling"] for t in fit.time)                  # CPU seconds, summed over chains
    return dict(mean=fit.stan_variable("log_T").mean(0), sd=fit.stan_variable("log_T").std(0),
                step=dr["stepsize__"].mean(), leap=dr["n_leapfrog__"].mean(), div=int(dr["divergent__"].sum()),
                tree=float((dr["treedepth__"] >= 10).mean()), ess=ess, rhat=float(cells["R_hat"].max()),
                per_grad=1000 * ess / grads, ms=1000 * t_samp / grads, per_sec=ess / t_samp, t_warm=t_warm, t_samp=t_samp)


def cols(r):
    return (f"{r['step']:.3g} | {r['leap']:.0f} | {r['div']} | {r['tree']:.0%} | {r['ess']:.0f} | {r['rhat']:.3f} | "
            f"{r['per_grad']:.2f} | {r['ms']:.3f} | {r['per_sec']:.1f}")


HEAD = "Step size | Leapfrogs | Divergences | At max depth | Min ESS | Max Rhat | ESS per 1,000 gradients | ms per gradient | ESS per second"


def main():
    rng = np.random.default_rng(20261005)
    T7 = load_table()
    T5 = collapse(T7, [[5], [1], [0], [2], [3, 4, 6]])             # SNP, Labour, Conservative, Lib Dem, others
    T3 = collapse(T7, [[5], [1], [0, 2, 3, 4, 6]])                 # SNP, Labour, others
    sizes = (("3x3", T3), ("5x5", T5), ("7x7", T7))
    L = ["# The single-area model in six coordinate systems", "",
         "One table, placeholder prior (ILR coordinates ~ normal(0, 2)), short runs. Times are CPU seconds in the "
         "sampling phase, summed over the two chains, on a 2-core cloud machine; treat them as relative.", "",
         "## Part 1: the stated Jacobian against finite differences", "",
         "Six random points per row. 'Offset' is the mean of (finite difference - stated); it should be zero, or a "
         "constant for the ILR parameterisation, whose Jacobian is stated up to a constant. 'Spread' is the largest "
         "deviation from that mean and should be near zero.", "",
         "| Size | Margin scaling | Parameterisation | Offset | Spread |", "|---|---|---|---|---|"]
    for label, T in (("3x3", T3), ("7x7", T7)):
        for scale in (0, 1):
            for name, off, spread in jacobian_check(T, 0.3, scale, rng):
                L.append(f"| {label} | {'on' if scale else 'off'} | {name} | {off:.2e} | {spread:.1e} |")
                print(L[-1], flush=True)

    L += ["", "## Part 2: sampling the same posterior", "",
          "2 chains, 500 warm-up + 500 draws each, started at the independence table, margin scaling off. "
          "Summaries are of the log cells; ESS is the smallest over cells, of 1,000 draws. 'Gap' is the largest "
          "difference in a posterior mean from the adjusted-row run, in posterior sds; 'sd ratio' is the range of "
          "the ratio of posterior sds to that run.", "",
          f"| Size | eps | Parameterisation | {HEAD} | Gap | sd ratio |", "|---|---|---|" + "---|" * 11]
    cost = {}
    for label, T in sizes:
        for eps in (1.0, 0.1):
            res = {p: sample(T, eps, p) for p, _ in PARAMS}
            for p, name in PARAMS:
                r = res[p]
                gap = (np.abs(r["mean"] - res[3]["mean"]) / res[3]["sd"]).max()
                ratio = r["sd"] / res[3]["sd"]
                L.append(f"| {label} | {eps:g} | {name} | {cols(r)} | {gap:.2f} | {ratio.min():.2f}-{ratio.max():.2f} |")
                print(L[-1], flush=True)
                cost.setdefault(p, {}).setdefault(label, []).append(r["ms"])

    L += ["", "## Part 3: longer runs of the sequential parameterisations at 7x7", "",
          "2 chains, 1,000 warm-up + 3,000 draws each, eps = 1, margin scaling off. ESS is of 6,000 draws; gap and "
          "sd ratio are against the adjusted-table run.", "",
          f"| Parameterisation | {HEAD} | Gap | sd ratio |", "|---|" + "---|" * 11]
    res = {p: sample(T7, 1.0, p, warm=1000, draws=3000) for p in (1, 2, 3, 4)}
    for p in (1, 2, 3, 4):
        r = res[p]
        gap = (np.abs(r["mean"] - res[4]["mean"]) / res[4]["sd"]).max()
        ratio = r["sd"] / res[4]["sd"]
        L.append(f"| {NAME[p]} | {cols(r)} | {gap:.2f} | {ratio.min():.2f}-{ratio.max():.2f} |")
        print(L[-1], flush=True)

    L += ["", "## Part 4: does scaling the margin parameters by the penalty width matter?", "",
          "Adjusted row, 2 chains, 500 warm-up + 500 draws. 'off': the margin parameters are plain log deviations "
          "from the observed margins. 'on': they are divided by the penalty's width, so they are about unit scale.", "",
          f"| Size | eps | Margin scaling | {HEAD} | Warm-up seconds |", "|---|---|---|" + "---|" * 10]
    for label, T in (("3x3", T3), ("7x7", T7)):
        for eps in (1.0, 0.1, 0.01):
            for scale in (0, 1):
                r = sample(T, eps, 3, scale=scale)
                L.append(f"| {label} | {eps:g} | {'on' if scale else 'off'} | {cols(r)} | {r['t_warm']:.1f} |")
                print(L[-1], flush=True)

    L += ["", "## Part 5: cost of one gradient by table size", "",
          "Milliseconds per gradient, mean of the two Part 2 runs. The last column is the 7x7 cost over the 3x3 cost; "
          "the number of parameters grows 5.4 times.", "",
          "| Parameterisation | 3x3 | 5x5 | 7x7 | 7x7 / 3x3 |", "|---|---|---|---|---|"]
    for p, name in PARAMS:
        v = [float(np.mean(cost[p][s])) for s in ("3x3", "5x5", "7x7")]
        L.append(f"| {name} | {v[0]:.3f} | {v[1]:.3f} | {v[2]:.3f} | {v[2] / v[0]:.1f} |")
        print(L[-1], flush=True)
    with open(os.path.join(ROOT, "experiments", "output", "02_single_area_models.md"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

import sys, numpy as np, pandas as pd
sys.path.insert(0, "experiments")
import importlib.util
spec = importlib.util.spec_from_file_location("r4", "experiments/04_real_table.py"); r4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r4)
e2 = r4.e2
T = r4.load("data/scot_2007_all_areas.csv"); R, C = T.shape; D = R*C
rows = []
for cen in ["indep", "diag", "reverse"]:
  for sb in [1.0, 2.0]:
    for p in [1, 2]:
      for seed in [1, 2]:
        d = e2.base_data(T, 0.1, 1)
        d.update(param=p, mu_b=d["V"].T @ r4.centre_table(T, cen).ravel(), sigma_b=np.full(D-1, sb), mu_logv=float(np.log(T.sum())), sigma_logv=1.0)
        f = e2.model.sample(data=d, chains=2, parallel_chains=2, iter_warmup=500, iter_sampling=500, seed=seed,
                            inits=[{"theta": e2.start(p, d)}]*2, show_progress=False, show_console=False)
        lt = f.stan_variable("log_T"); Tt = np.exp(lt)
        w = Tt.sum(2); m = Tt.sum(1); tot = Tt.sum((1, 2))
        lo = np.maximum(0, w[:, :, None] + m[:, None, :] - tot[:, None, None]); hi = np.minimum(w[:, :, None], m[:, None, :])
        pos = (Tt - lo) / (hi - lo)
        s = f.summary(); c = s[s.index.str.startswith("log_T[")]
        ess = c["ESS_bulk"].values.reshape(R, C) if "ESS_bulk" in c else c["N_Eff"].values.reshape(R, C)
        mp = pos.mean(0); wd = (hi - lo).mean(0) / hi.mean(0)
        near = ((mp < 0.02) | (mp > 0.98)).sum()
        rows.append(dict(centre=cen, sigma_b=sb, param=r4.NAMES[p], seed=seed, rhat=c["R_hat"].max(), min_ess=ess.min(),
                         cells_near_bound=int(near), frac_draws_near=float(((pos < 0.02) | (pos > 0.98)).mean()),
                         pos_sd_median=float(np.median(pos.std(0))), corr_ess_minpos=float(np.corrcoef(ess.ravel(), np.minimum(mp, 1-mp).ravel())[0,1]),
                         corr_ess_width=float(np.corrcoef(ess.ravel(), wd.ravel())[0,1]), worst_cell_width=float(wd.ravel()[ess.argmin()])))
        print(rows[-1], flush=True)
df = pd.DataFrame(rows); df.to_csv("experiments/output/05_bounds.csv", index=False)
print(df.groupby(["centre","sigma_b","param"]).agg(rhat=("rhat","max"), min_ess=("min_ess","median"), near=("cells_near_bound","median"), frac=("frac_draws_near","median"), posd=("pos_sd_median","median"), cmin=("corr_ess_minpos","median"), cw=("corr_ess_width","median")).round(3).to_string())

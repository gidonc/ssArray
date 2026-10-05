# Stage A evaluation in R: tables simulated from the model's own prior.
# Mirrors experiments/03_simulated_tables.py (which has been run); THIS FILE HAS NOT BEEN RUN YET.
#
# For each simulated table: ILR coordinates ~ normal(0, sigma_b), log total ~ normal(log total, 1);
# the row and column totals are observed with the model's own noise, obs = true + eps * sqrt(true) * N(0, 1);
# the same data are then fitted in each parameterisation. Because the data come from the model, the rank of
# the true log cell among the posterior draws should be uniform (simulation-based calibration).
#
# Needs R/single_area.R sourced first (sa_helmert, sa_data, sa_init).

PARAMS <- c("ILR + log volume" = 0, "log cells" = 5, "position logit" = 1,
            "odds-ratio logit" = 2, "adjusted row" = 3, "adjusted table" = 4)

# one row per task: n_rep simulated tables at each (sigma_b, eps); a task fits all parameterisations to one table
sim_design <- function(n_rep = 30, sigma_b = c(1, 2, 4), eps = c(1, 0.1, 0.01)) {
  g <- expand.grid(rep = seq_len(n_rep), eps = eps, sigma_b = sigma_b)
  g$task <- seq_len(nrow(g))
  g[, c("task", "sigma_b", "eps", "rep")]
}

sim_table <- function(sigma_b, eps, R = 7, C = 7, total = 20000) {
  D  <- R * C
  V  <- sa_helmert(D)
  z  <- rnorm(D - 1, 0, sigma_b)
  lt <- as.vector(V %*% z)
  lt <- lt + rnorm(1, log(total), 1) - log(sum(exp(lt)))
  T  <- matrix(exp(lt), R, C, byrow = TRUE)            # cells were ordered row by row
  w  <- rowSums(T); m <- colSums(T)
  list(T = T,
       obs_w = pmax(w + eps * sqrt(w) * rnorm(R), 1),
       obs_m = pmax(m + eps * sqrt(m) * rnorm(C), 1))
}

# fit one parameterisation to one simulated table; returns list(row = one-row data.frame, ranks = vector)
sim_fit_one <- function(mod, sim, param, eps, sigma_b, seed, total = 20000,
                        chains = 2, warmup = 500, draws = 500) {
  d <- sa_data(sim$T, param = param, eps = eps, sigma_b = sigma_b)
  d$obs_w     <- sim$obs_w
  d$obs_m     <- sim$obs_m
  d$row_order <- order(d$obs_w)
  d$mu_logv   <- log(total)
  t0  <- Sys.time()
  fit <- mod$sample(data = d, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = draws,
                    init = rep(list(sa_init(d)), chains), seed = seed,
                    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  s    <- fit$summary("log_T")
  dg   <- fit$diagnostic_summary(quiet = TRUE)
  sd_  <- fit$sampler_diagnostics(format = "df")
  lT   <- as.matrix(fit$draws("log_T", format = "draws_matrix"))   # columns column-major, like as.vector(T)
  ranks <- colMeans(sweep(lT, 2, as.vector(log(sim$T)), "<"))
  list(row = data.frame(min_ess = min(s$ess_bulk), max_rhat = max(s$rhat),
                        divergences = sum(dg$num_divergent), max_treedepth = sum(dg$num_max_treedepth),
                        min_ebfmi = min(dg$ebfmi), leapfrogs = sum(sd_$n_leapfrog__),
                        seconds = wall, rank_mean = mean(ranks), rank_sd = sd(ranks), error = NA_character_),
       ranks = ranks)
}

sim_fit_safe <- function(...) {
  tryCatch(sim_fit_one(...), error = function(e)
    list(row = data.frame(min_ess = NA_real_, max_rhat = NA_real_, divergences = NA_integer_, max_treedepth = NA_integer_,
                          min_ebfmi = NA_real_, leapfrogs = NA_real_, seconds = NA_real_, rank_mean = NA_real_,
                          rank_sd = NA_real_, error = substr(conditionMessage(e), 1, 200)),
         ranks = rep(NA_real_, 49)))
}

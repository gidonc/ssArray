# Run the single-area model (stan/single_area.stan) from R with cmdstanr.
#
# NOT YET RUN: written without access to R. It mirrors experiments/02_single_area_models.py, which has been run.
#
# Usage, from the repository root:
#   source("R/single_area.R")
#   tab <- read.csv("data/scot_area60.csv")
#   T7  <- xtabs(votes ~ row_no + col_no, tab)            # 7 x 7 table of the example area
#   fit <- sa_fit(T7, param = 3, eps = 1)                 # 3 = adjusted row
#   fit$summary("log_T")
#
# param: 0 = ILR + log volume, 5 = log cells,
#        1 = between bounds, 2 = cell logit, 3 = adjusted row, 4 = adjusted table
# eps:   relative width of the margin penalty (1 is about Poisson noise; smaller is tighter)

# orthonormal basis of the contrasts among D cells (Helmert)
sa_helmert <- function(D) {
  V <- matrix(0, D, D - 1)
  for (k in 1:(D - 1)) {
    V[1:k, k] <- 1 / sqrt(k * (k + 1))
    V[k + 1, k] <- -k / sqrt(k * (k + 1))
  }
  V
}

# data list for the model; only the margins of `tab` are used
sa_data <- function(tab, param, eps, sigma_b = 2, delta = 0) {
  tab <- unclass(as.matrix(tab))
  R <- nrow(tab); C <- ncol(tab); D <- R * C
  w <- as.numeric(rowSums(tab)); m <- as.numeric(colSums(tab))
  list(
    R = R, C = C, param = param, obs_w = w, obs_m = m, eps = eps,
    V = sa_helmert(D),
    mu_b = rep(0, D - 1), sigma_b = rep(sigma_b, D - 1),      # placeholder prior
    mu_logv = log(sum(w)), sigma_logv = 1,
    row_order = order(w),                                     # smallest row first, largest found by subtraction
    rem_col = pmin(seq_len(R), C),                            # each row's own column is its reference
    delta = delta
  )
}

# starting values on the margins: the independence table
sa_init <- function(d) {
  T0 <- outer(d$obs_w, d$obs_m) / sum(d$obs_w)
  lt <- log(as.vector(t(T0)))                                 # cells row by row
  theta <- if (d$param == 0) c(as.vector(t(d$V) %*% lt), log(sum(T0)))
           else if (d$param == 5) lt
           else rep(0, d$R * d$C)
  list(theta = theta)
}

sa_fit <- function(tab, param, eps, chains = 2, iter_warmup = 500, iter_sampling = 500, ...) {
  d <- sa_data(tab, param, eps)
  mod <- cmdstanr::cmdstan_model("stan/single_area.stan", include_paths = "stan")
  mod$sample(data = d, chains = chains, parallel_chains = chains,
             iter_warmup = iter_warmup, iter_sampling = iter_sampling,
             init = rep(list(sa_init(d)), chains), ...)
}

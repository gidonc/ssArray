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

# Margin coordinates for the sequential versions.
# The model's margin parameters are z = (z_row[R], z_col[C-1]): row totals obs_w * exp(s_r z_row), column shares a softmax
# of (observed log-ratios against the LAST column + s_c z_col). The data matrix K_margin sets z = K_margin %*% u, and u is
# what is sampled. Any K is a fixed linear change of coordinates, so the density on tables is the same.
#   "last"      identity: the original coordinates
#   "largest"   column log-ratios against the largest column instead of the last
#   "whitened"  u roughly unit scale and uncorrelated under the margin penalty, no reference column (uses scale_margins = 0)
# A small reference column puts its noise into every column parameter, which a diagonal metric cannot undo
# (experiments/21_reference_column.py), so "whitened" or "largest" is the safer choice when the last column is small.
# Returns list(K = matrix, scale_margins = value to pass to the model).
sa_margin_K <- function(w, m, eps, kind = c("last", "largest", "whitened"), scale_margins = 1) {
  kind <- match.arg(kind)
  R <- length(w); C <- length(m); n <- R + C - 1
  if (kind == "last") return(list(K = diag(n), scale_margins = scale_margins))
  if (kind == "largest") {
    j <- which.max(m)
    if (j == C) return(list(K = diag(n), scale_margins = scale_margins))
    s_old <- if (scale_margins == 1) eps * sqrt(1 / m[-C] + 1 / m[C]) else rep(1, C - 1)
    others <- setdiff(seq_len(C), j)                          # new parameters: columns other than j, in natural order
    t_new <- if (scale_margins == 1) eps * sqrt(1 / m[others] + 1 / m[j]) else rep(1, C - 1)
    # e_c = log-ratio deviation of column c against column j; old d_c (against the last) = e_c - e_last, d_j = -e_last
    A <- matrix(0, C - 1, C - 1); last <- match(C, others)
    for (c in seq_len(C - 1)) {
      if (c != j) A[c, match(c, others)] <- A[c, match(c, others)] + 1
      A[c, last] <- A[c, last] - 1
    }
    K <- diag(n)
    K[(R + 1):n, (R + 1):n] <- sweep(A, 2, t_new, "*") / s_old
    return(list(K = K, scale_margins = scale_margins))
  }
  N <- sum(w); p <- w / N; q <- m / N
  M <- matrix(0, R + C, n)                                    # d (log w, log m) / d z, unit scales
  M[1:R, 1:R] <- diag(R)
  M[(R + 1):(R + C), 1:R] <- matrix(p, C, R, byrow = TRUE)    # the total moves every column
  B <- rbind(diag(C - 1), 0) - matrix(q[-C], C, C - 1, byrow = TRUE)
  M[(R + 1):(R + C), (R + 1):n] <- B
  H <- t(M) %*% diag(c(w, m) / eps^2) %*% M + diag(n)         # penalty precision, + 1 so no parameter is wider than sd 1
  L <- t(chol(H))                                             # H = L L'
  list(K = t(solve(L)), scale_margins = 0)
}

# data list for the model; only the margins of `tab` are used
sa_data <- function(tab, param, eps, sigma_b = 2, delta = 0, scale_margins = 0,
                    margin_coords = c("last", "largest", "whitened")) {
  tab <- unclass(as.matrix(tab))
  R <- nrow(tab); C <- ncol(tab); D <- R * C
  w <- as.numeric(rowSums(tab)); m <- as.numeric(colSums(tab))
  mk <- sa_margin_K(w, m, eps, match.arg(margin_coords), scale_margins)
  list(
    R = R, C = C, param = param, obs_w = w, obs_m = m, eps = eps,
    scale_margins = mk$scale_margins,                         # sequential versions: 1 = margin parameters on the penalty's scale
    K_margin = mk$K,                                          # sequential versions: linear change of margin coordinates
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

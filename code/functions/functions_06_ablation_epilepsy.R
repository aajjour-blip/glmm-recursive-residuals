# ============================================================
# Function definitions for code/06_ablation_epilepsy.R
# Sourced by that script; run it from the root of the supplement.
# ============================================================

# ----------------------------------------------------------------------
# fit_model(): Fit the glmmTMB model described by a specification to a data frame.
# Arguments:
#   spec: model specification list (label, family, formula, fixed, index, seeds)
#   d: data frame of residuals and index variables
# Returns: fitted glmmTMB object
fit_model <- function(spec, d) {
  fam <- if (spec$family == "poisson") stats::poisson() else glmmTMB::nbinom2()
  suppressWarnings(glmmTMB::glmmTMB(spec$formula, data = d, family = fam,
                                    control = glmmTMB::glmmTMBControl(parallel = 1)))
}

# ----------------------------------------------------------------------
# re_sd(): Random-intercept standard deviation of a glmmTMB fit, bounded below by 1e-8.
# Arguments:
#   fit: fitted GLMM object
# Returns: numeric scalar
re_sd <- function(fit) {
  s <- as.numeric(attr(VarCorr(fit)$cond$subject, "stddev"))[1]
  if (!is.finite(s) || s < 1e-8) s <- 1e-8
  s
}

# ----------------------------------------------------------------------
# nb2_theta(): NB2 size parameter of a glmmTMB fit (sigma(fit)).
# Arguments:
#   fit: fitted GLMM object
# Returns: numeric scalar
nb2_theta <- function(fit) {
  th <- as.numeric(sigma(fit))
  if (!is.finite(th) || th <= 0) stop("Could not extract a valid NB2 dispersion via sigma(fit).")
  th
}

# ----------------------------------------------------------------------
# cdf_count(): Poisson or NB2 distribution function.
# Arguments:
#   k: observed count(s)
#   lam: conditional mean(s)
#   family: "poisson" or "nbinom2"
#   theta: NB2 size (dispersion) parameter; ignored for Poisson
# Returns: numeric vector of probabilities
cdf_count <- function(k, lam, family, theta) {
  if (family == "poisson") stats::ppois(k, lam) else stats::pnbinom(k, mu = lam, size = theta)
}

# ----------------------------------------------------------------------
# pmf_count(): Poisson or NB2 probability mass function.
# Arguments:
#   k: observed count(s)
#   lam: conditional mean(s)
#   family: "poisson" or "nbinom2"
#   theta: NB2 size (dispersion) parameter; ignored for Poisson
# Returns: numeric vector of probabilities
pmf_count <- function(k, lam, family, theta) {
  if (family == "poisson") stats::dpois(k, lam) else stats::dnbinom(k, mu = lam, size = theta)
}

# ----------------------------------------------------------------------
# safe_qnorm(): Standard normal quantile function with probabilities truncated to [eps, 1 -
#  eps].
# Arguments:
#   u: numeric vector of probabilities
# Returns: numeric vector of normal quantiles
safe_qnorm <- function(u) { eps <- .Machine$double.eps; stats::qnorm(pmin(pmax(u, eps), 1 - eps)) }

# ----------------------------------------------------------------------
# colcumsum(): Column-wise cumulative sums of a matrix.
# Arguments:
#   M: numeric matrix
# Returns: matrix of the same dimensions
colcumsum <- function(M) {
  if (nrow(M) > 1) for (t in 2:nrow(M)) M[t, ] <- M[t - 1, ] + M[t, ]
  M
}

# ----------------------------------------------------------------------
# rr_constructions(): Recursive residuals under the sequential, marginal and full-cluster
#  posterior constructions (Section 6.6), using common uniform draws.
# Arguments:
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   spec: model specification list (label, family, formula, fixed, index, seeds)
#   seed: random-number seed
#   which: residual constructions to compute (subset of "seq", "marg", "full")
# Returns: data frame sorted by subject and period with columns rr_<w> and mean_<w> for each
#  construction
# seq : weights proportional to prior x likelihood of y_i1..y_i,t-1   (proposed, Eq. 5)
# marg: weights equal to the prior of b_i                            (no conditioning)
# full: weights proportional to prior x likelihood of y_i1..y_i,n_i  (full-cluster posterior)
# All constructions use the same uniform draws a_it, drawn in the same
# order as in scripts 03 and 05, so the "seq" residuals are identical to
# the published ones for the same seed.
rr_constructions <- function(fit, data, spec, seed, which = CONSTR) {
  set.seed(seed)
  d <- data[order(data$subject, data$period_num), ]
  rownames(d) <- NULL
  beta <- fixef(fit)$cond
  X <- stats::model.matrix(spec$fixed, d)
  miss <- setdiff(names(beta), colnames(X))
  if (length(miss) > 0) stop("Missing fixed-effect columns: ", paste(miss, collapse = ", "))
  eta <- as.numeric(X[, names(beta), drop = FALSE] %*% beta)
  theta <- if (spec$family == "nbinom2") nb2_theta(fit) else NA_real_
  gh <- statmod::gauss.quad.prob(n = N_QUAD, dist = "normal")
  bq <- re_sd(fit) * gh$nodes
  lw0 <- log(gh$weights)
  Q <- length(bq)
  n <- nrow(d)
  a <- stats::runif(n)
  y <- d$y
  RR <- matrix(NA_real_, n, length(which), dimnames = list(NULL, which))
  MU <- RR
  ids <- d$subject
  for (id in unique(ids)) {
    rows <- which(ids == id)
    ni <- length(rows)
    lam <- exp(outer(eta[rows], bq, "+"))
    kk <- rep(y[rows], times = Q)
    lv <- as.vector(lam)
    Glo <- matrix(cdf_count(kk - 1, lv, spec$family, theta), nrow = ni)
    Ghi <- matrix(cdf_count(kk, lv, spec$family, theta), nrow = ni)
    L <- log(pmax(matrix(pmf_count(kk, lv, spec$family, theta), nrow = ni), .Machine$double.xmin))
    prev <- rbind(0, colcumsum(L))[seq_len(ni), , drop = FALSE]
    tot <- colSums(L)
    for (t in seq_len(ni)) {
      row <- rows[t]
      for (w in which) {
        lw <- switch(w, seq = lw0 + prev[t, ], marg = lw0, full = lw0 + tot)
        pw <- exp(lw - max(lw)); pw <- pw / sum(pw)
        glo <- sum(pw * Glo[t, ]); ghi <- sum(pw * Ghi[t, ])
        RR[row, w] <- safe_qnorm(glo + a[row] * (ghi - glo))
        MU[row, w] <- sum(pw * lam[t, ])
      }
    }
  }
  for (w in which) { d[[paste0("rr_", w)]] <- RR[, w]; d[[paste0("mean_", w)]] <- MU[, w] }
  d
}

# ----------------------------------------------------------------------
# index_values(): Values of an index for a construction (the fitted-mean index depends on
#  the construction).
# Arguments:
#   d: data frame of residuals and index variables
#   w: residual construction key ("seq", "marg" or "full")
#   v: index name, or numeric vector of reference values (sim_pvalues)
# Returns: vector of index values
index_values <- function(d, w, v) if (v == "pred_mean") d[[paste0("mean_", w)]] else d[[v]]

# ----------------------------------------------------------------------
# proc_W(): Cumulative residual process for one construction and one index.
# Arguments:
#   d: data frame of residuals and index variables
#   w: residual construction key ("seq", "marg" or "full")
#   v: index name, or numeric vector of reference values (sim_pvalues)
# Returns: numeric vector of process values
proc_W <- function(d, w, v) {
  q <- index_values(d, w, v)
  r <- d[[paste0("rr_", w)]]
  o <- order(q, seq_len(nrow(d)))
  cumsum(r[o]) / sqrt(length(r))
}

# ----------------------------------------------------------------------
# stats_one(): Component statistics for one construction, plus the direct statistics T_KS
#  and T_CvM.
# Arguments:
#   d: data frame of residuals and index variables
#   w: residual construction key ("seq", "marg" or "full")
#   index_vars: character vector of index columns for the cumulative processes
# Returns: named numeric vector of statistics
stats_one <- function(d, w, index_vars) {
  r <- d[[paste0("rr_", w)]]
  N <- length(r)
  out <- numeric(0)
  for (v in index_vars) {
    if (length(unique(index_values(d, w, v))) < 2) next
    W <- proc_W(d, w, v)
    out[paste0("Tinf_", v)] <- max(abs(W))
    out[paste0("T2_", v)] <- mean(W^2)
  }
  out["Tmax_r"] <- max(abs(r))
  C <- tapply(r, d$subject, function(z) sum(z) / sqrt(length(z)))
  out["Tmax_C"] <- max(abs(C), na.rm = TRUE)
  out["Tdisp"] <- abs(mean(r^2) - 1)
  same <- c(FALSE, as.character(d$subject[-1]) == as.character(d$subject[-N]))
  out["Tlag"] <- if (sum(same) > 3) abs(stats::cor(r[same], r[which(same) - 1L])) else NA_real_
  # direct distributional statistics on the same residuals
  rs <- sort(r); Fh <- stats::pnorm(rs); i <- seq_len(N)
  out["Tks"] <- max(pmax(i / N - Fh, Fh - (i - 1) / N))
  out["Tcvm"] <- 1 / (12 * N) + sum((Fh - (2 * i - 1) / (2 * N))^2)
  out
}

# ----------------------------------------------------------------------
# p_upper(): Upper-tail calibrated p-value of an observed statistic.
# Arguments:
#   obs: observed value of a statistic
#   sim: numeric vector of reference values of the same statistic
# Returns: p-value, or NA
p_upper <- function(obs, sim) {
  s <- sim[is.finite(sim)]
  if (!is.finite(obs) || length(s) == 0) return(NA_real_)
  (1 + sum(s >= obs)) / (length(s) + 1)
}

# ----------------------------------------------------------------------
# sim_pvalues(): Calibrated p-value of each reference value within its own reference sample.
# Arguments:
#   v: index name, or numeric vector of reference values (sim_pvalues)
# Returns: numeric vector of p-values
sim_pvalues <- function(v) {
  p <- rep(NA_real_, length(v)); ok <- is.finite(v); n <- sum(ok)
  if (n > 0) { rk <- rank(v[ok], ties.method = "min"); p[ok] <- (n + 2 - rk) / (n + 1) }
  p
}

# ----------------------------------------------------------------------
# calib(): Component p-values and calibrated minimum-p omnibus p-value over a set of
#  statistics.
# Arguments:
#   obs: observed value of a statistic
#   sim: numeric vector of reference values of the same statistic
#   set: names of the statistics combined in the omnibus test
# Returns: list with p, p_min, p_combo, lead and B
calib <- function(obs, sim, set) {
  set <- set[set %in% names(obs) & set %in% colnames(sim)]
  set <- set[is.finite(obs[set])]
  p <- vapply(set, function(s) p_upper(obs[[s]], sim[, s]), numeric(1))
  p_min <- min(p, na.rm = TRUE)
  sp <- matrix(vapply(set, function(s) sim_pvalues(sim[, s]), numeric(nrow(sim))), nrow = nrow(sim))
  sim_min <- suppressWarnings(apply(sp, 1, min, na.rm = TRUE))
  list(p = p, p_min = p_min,
       p_combo = (1 + sum(sim_min <= p_min)) / (length(sim_min) + 1),
       lead = if (any(is.finite(p))) set[which.min(p)] else NA_character_, B = nrow(sim))
}

# ----------------------------------------------------------------------
# to_mat(): Bind a list of statistic vectors into a matrix with fixed column names.
# Arguments:
#   lst: list of named numeric vectors of statistics
#   nm: names of the statistics (columns of the result)
# Returns: numeric matrix, one row per replicate
to_mat <- function(lst, nm) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  m <- do.call(rbind, lapply(lst, function(v) unname(v[nm])))
  colnames(m) <- nm
  m
}

# ----------------------------------------------------------------------
# diagnose(): Run all ablation variants for one fitted model and dataset (Section 6.6,
#  Appendix B).
# Arguments:
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   spec: model specification list (label, family, formula, fixed, index, seeds)
#   B: number of simulated reference datasets
#   seeds: list of seeds: obs, sim, res and kp
#   keep_process: index names whose simulated processes are stored for plotting
# Returns: list with variant results V, classical p-values, observed statistics, residuals,
#  stored processes and failure counts
diagnose <- function(fit, data, spec, B, seeds, keep_process = character(0)) {
  d_obs <- rr_constructions(fit, data, spec, seeds$obs)
  obs <- lapply(stats::setNames(CONSTR, CONSTR), function(w) stats_one(d_obs, w, spec$index))
  classical <- c(ks = suppressWarnings(stats::ks.test(d_obs$rr_seq, "pnorm")$p.value),
                 sw = stats::shapiro.test(d_obs$rr_seq)$p.value)
  S <- list(seq = vector("list", B), marg = vector("list", B), full = vector("list", B),
            plugin = vector("list", B), plugin_marg = vector("list", B), known = vector("list", B))
  P <- list()
  for (w in CONSTR) for (v in keep_process) P[[paste(w, v, sep = ":")]] <- matrix(NA_real_, B,
      nrow(d_obs))
  n_refit_fail <- 0L; n_skip <- 0L
  for (b in seq_len(B)) {
    ysim <- tryCatch(as.numeric(stats::simulate(fit, nsim = 1, seed = seeds$sim + b)[[1]]),
                     error = function(e) NULL)
    if (is.null(ysim)) { n_skip <- n_skip + 1L; next }
    dat_b <- data; dat_b$y <- ysim
    fit_b <- tryCatch(fit_model(spec, dat_b), error = function(e) NULL)
    if (is.null(fit_b)) { n_refit_fail <- n_refit_fail + 1L; fit_b <- fit }
    d_b <- tryCatch(rr_constructions(fit_b, dat_b, spec, seeds$res + b), error = function(e) NULL)
    d_p <- tryCatch(rr_constructions(fit, dat_b, spec, seeds$res + b, which = c("seq",
        "marg")), error = function(e) NULL)
    if (is.null(d_b) || is.null(d_p)) { n_skip <- n_skip + 1L; next }
    for (w in CONSTR) S[[w]][[b]] <- stats_one(d_b, w, spec$index)
    S$plugin[[b]] <- stats_one(d_p, "seq", spec$index)
    S$plugin_marg[[b]] <- stats_one(d_p, "marg", spec$index)
    set.seed(seeds$kp + b)
    d_k <- d_obs
    d_k$rr_seq <- stats::rnorm(nrow(d_k))
    S$known[[b]] <- stats_one(d_k, "seq", spec$index)
    for (w in CONSTR) for (v in keep_process) P[[paste(w, v, sep = ":")]][b, ] <- proc_W(d_b, w, v)
  }
  nm <- names(obs$seq)
  M <- lapply(S, to_mat, nm = nm)
  idx_names <- grep("^T(inf|2)_", nm, value = TRUE)
  SETS <- list(full = c(idx_names, "Tmax_r", "Tmax_C", "Tdisp", "Tlag"),
               indexed = idx_names,
               direct = c("Tks", "Tcvm", "Tdisp", "Tmax_r"))
  V <- list(proposed = calib(obs$seq,  M$seq,    SETS$full),
            marginal = calib(obs$marg, M$marg,   SETS$full),
            fullpost = calib(obs$full, M$full,   SETS$full),
            plugin   = calib(obs$seq,  M$plugin, SETS$full),
            known    = calib(obs$seq,  M$known,  SETS$full),
            neither  = calib(obs$marg, M$plugin_marg, SETS$full),
            indexed  = calib(obs$seq,  M$seq,    SETS$indexed),
            direct   = calib(obs$seq,  M$seq,    SETS$direct))
  list(V = V, classical = classical, obs = obs, d_obs = d_obs, P = P,
       n_refit_fail = n_refit_fail, n_skip = n_skip)
}

# ----------------------------------------------------------------------
# summ_row(): Collect the omnibus p-value and leading component of every variant.
# Arguments:
#   res: list returned by diagnose()
# Returns: data frame with columns key, p and lead
summ_row <- function(res) {
  p <- c(vapply(res$V, function(x) x$p_combo, numeric(1)),
         ks = unname(res$classical["ks"]), sw = unname(res$classical["sw"]))
  lead <- c(vapply(res$V, function(x) x$lead, character(1)), ks = NA_character_, sw = NA_character_)
  data.frame(key = VARIANTS$key, p = unname(p[VARIANTS$key]), lead = unname(lead[VARIANTS$key]),
             stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------------
# stat_tex(): LaTeX label of a statistic name.
# Arguments:
#   s: statistic name
# Returns: character string
stat_tex <- function(s) {
  fixed_map <- c(Tmax_r = "$T_{\\max r}$", Tmax_C = "$T_{\\max C}$", Tdisp = "$T_{\\mathrm{disp}}$",
                 Tlag = "$T_{\\mathrm{lag}}$", Tks = "$T_{\\mathrm{KS}}$",
                     Tcvm = "$T_{\\mathrm{CvM}}$")
  if (s %in% names(fixed_map)) return(fixed_map[[s]])
  idx_lab <- c(period_num = "period", pred_mean = "fitted mean", base = "baseline",
               log_base = "log baseline", age = "age", trt = "treatment")
  m <- regmatches(s, regexec("^(Tinf|T2)_(.*)$", s))[[1]]
  if (length(m) == 3) {
    lab <- if (m[3] %in% names(idx_lab)) idx_lab[[m[3]]] else m[3]
    sym <- if (m[2] == "Tinf") "\\infty" else "2"
    return(paste0("$T_{", sym, "}$ (", lab, ")"))
  }
  s
}

# ----------------------------------------------------------------------
# fmt3(): Format p-values with three decimals.
# Arguments:
#   p: proportion(s)
# Returns: character vector
fmt3 <- function(p) ifelse(is.finite(p), sprintf("%.3f", p), "--")

# ----------------------------------------------------------------------
# block_col(): Blank repeated block labels for table output.
# Arguments:
#   blocks: character vector of block labels
# Returns: character vector
block_col <- function(blocks) {
  out <- blocks
  if (length(blocks) > 1) for (i in 2:length(blocks)) if (blocks[i] == blocks[i - 1]) out[i] <- ""
  out
}

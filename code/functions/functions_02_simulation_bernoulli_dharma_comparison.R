# ============================================================
# Function definitions for code/02_simulation_bernoulli_dharma_comparison.R
# Sourced by that script; run it from the root of the supplement.
# ============================================================

# ----------------------------------------------------------------------
# load_required_packages(): Stop with an informative message if any required package is not
#  installed.
# Arguments:
#   pkgs: character vector of package names
# Returns: invisible TRUE
load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop("Install the following R packages before running this script: ", paste(missing,
        collapse = ", "))
  }
  invisible(TRUE)
}

# ----------------------------------------------------------------------
# as_int_env(): Read an integer setting from an environment variable.
# Arguments:
#   name: name of the environment variable
#   default: value returned when the variable is unset or invalid
# Returns: integer value
as_int_env <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) return(as.integer(default))
  out <- suppressWarnings(as.integer(val))
  if (!is.finite(out)) as.integer(default) else out
}

# ----------------------------------------------------------------------
# as_num_env(): Read a numeric setting from an environment variable.
# Arguments:
#   name: name of the environment variable
#   default: value returned when the variable is unset or invalid
# Returns: numeric value
as_num_env <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(val))
  if (!is.finite(out)) as.numeric(default) else out
}

# ----------------------------------------------------------------------
# as_logical_env(): Read a logical setting from an environment variable.
# Arguments:
#   name: name of the environment variable
#   default: value returned when the variable is unset or invalid
# Returns: logical value
as_logical_env <- function(name, default = FALSE) {
  val <- toupper(Sys.getenv(name, unset = if (default) "TRUE" else "FALSE"))
  val %in% c("TRUE", "T", "1", "YES", "Y")
}

# ----------------------------------------------------------------------
# normalize_log_weights(): Convert log weights to normalized weights, stabilized by
#  subtracting the maximum.
# Arguments:
#   logw: numeric vector of log weights at the quadrature nodes
# Returns: numeric vector of weights summing to one
normalize_log_weights <- function(logw) {
  m <- max(logw)
  w <- exp(logw - m)
  w / sum(w)
}

# ----------------------------------------------------------------------
# safe_qnorm(): Standard normal quantile function with probabilities truncated to [eps, 1 -
#  eps].
# Arguments:
#   u: numeric vector of probabilities
# Returns: numeric vector of normal quantiles
safe_qnorm <- function(u) {
  eps <- .Machine$double.eps
  qnorm(pmin(pmax(u, eps), 1 - eps))
}

# ----------------------------------------------------------------------
# rr_process(): Cumulative residual process W_N(q) of Eq. (10), with observations ordered by
#  the index.
# Arguments:
#   rr_data: data frame of recursive residuals and index variables (one row per observation)
#   q_var: name of the column used as the index of the cumulative process
#   rr_var: name of the residual column (default "rr")
# Returns: data frame with columns q (sorted index) and W (process values)
rr_process <- function(rr_data, q_var, rr_var = "rr") {
  ord <- order(rr_data[[q_var]], seq_len(nrow(rr_data)))
  q_sorted <- rr_data[[q_var]][ord]
  r_sorted <- rr_data[[rr_var]][ord]
  W <- cumsum(r_sorted) / sqrt(length(r_sorted))
  data.frame(q = q_sorted, W = W)
}

# ----------------------------------------------------------------------
# component_stats(): Compute the component statistics of Section 4.2: T_inf and T_2 for each
#  index, T_max r, T_max C, T_disp and T_lag.
# Arguments:
#   rr_data: data frame of recursive residuals and index variables (one row per observation)
#   index_vars: character vector of index columns for the cumulative processes
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   rr_var: name of the residual column (default "rr")
# Returns: one-row data frame of component statistics
component_stats <- function(rr_data, index_vars, id_var, order_var, rr_var = "rr") {
  out <- list()
  for (v in index_vars) {
    if (!v %in% names(rr_data)) next
    if (all(is.na(rr_data[[v]]))) next
    if (length(unique(rr_data[[v]])) < 2) next
    proc <- rr_process(rr_data, v, rr_var)
    nm <- make.names(v)
    out[[paste0("Tinf_", nm)]] <- max(abs(proc$W), na.rm = TRUE)
    out[[paste0("T2_", nm)]] <- mean(proc$W^2, na.rm = TRUE)
  }

  r <- rr_data[[rr_var]]
  out[["Tmax_r"]] <- max(abs(r), na.rm = TRUE)

  C <- rr_data |>
    group_by(.data[[id_var]]) |>
    summarise(C = sum(.data[[rr_var]]) / sqrt(n()), .groups = "drop")

  out[["Tmax_C"]] <- max(abs(C$C), na.rm = TRUE)
  out[["Tdisp"]] <- abs(mean(r^2, na.rm = TRUE) - 1)

  rr_lag <- rr_data |>
    arrange(.data[[id_var]], .data[[order_var]]) |>
    group_by(.data[[id_var]]) |>
    mutate(rr_lag = lag(.data[[rr_var]])) |>
    ungroup() |>
    filter(!is.na(rr_lag))

  out[["Tlag"]] <- if (nrow(rr_lag) > 3) {
    abs(cor(rr_lag[[rr_var]], rr_lag$rr_lag))
  } else {
    NA_real_
  }

  as.data.frame(out)
}

# ----------------------------------------------------------------------
# p_from_sim(): Refit-calibrated upper-tail p-value of an observed statistic, (1 + #{sim >=
#  obs}) / (B + 1).
# Arguments:
#   obs: observed value of a statistic
#   sim: numeric vector of reference values of the same statistic
# Returns: p-value, or NA if no finite reference values
p_from_sim <- function(obs, sim) {
  valid <- sim[is.finite(sim)]
  if (!is.finite(obs) || length(valid) == 0) return(NA_real_)
  (1 + sum(valid >= obs, na.rm = TRUE)) / (length(valid) + 1)
}

# ----------------------------------------------------------------------
# calibrate_components(): Component p-values and the calibrated minimum-p omnibus p-value
#  p_combo of Eq. (25).
# Arguments:
#   obs_stats: one-row data frame of observed component statistics
#   sim_stats: data frame of component statistics, one row per reference replicate
# Returns: data frame with component p-values, p_min, p_combo, leading statistic and
#  detection flag
calibrate_components <- function(obs_stats, sim_stats) {
  stat_names <- intersect(names(obs_stats), names(sim_stats))
  stat_names <- stat_names[sapply(obs_stats[stat_names], function(x) is.finite(as.numeric(x[1])))]
  stat_names <- stat_names[sapply(sim_stats[stat_names], function(x) any(is.finite(x)))]

  if (length(stat_names) == 0) {
    return(data.frame(p_combo_minP = NA_real_, p_min_raw = NA_real_, best_stat = NA_character_))
  }

  raw_p <- sapply(stat_names, function(s) p_from_sim(obs_stats[[s]][1], sim_stats[[s]]))
  p_min_raw <- min(raw_p, na.rm = TRUE)

  sim_p <- matrix(NA_real_, nrow = nrow(sim_stats), ncol = length(stat_names))
  colnames(sim_p) <- stat_names

  for (j in seq_along(stat_names)) {
    s <- stat_names[j]
    vals <- sim_stats[[s]]
    sim_p[, j] <- sapply(vals, function(v) p_from_sim(v, vals))
  }

  sim_min <- apply(sim_p, 1, min, na.rm = TRUE)
  p_combo <- (1 + sum(sim_min <= p_min_raw, na.rm = TRUE)) / (length(sim_min) + 1)

  data.frame(
    p_combo_minP = p_combo,
    p_min_raw = p_min_raw,
    best_stat = names(which.min(raw_p))
  )
}

# ----------------------------------------------------------------------
# mc_se(): Monte Carlo standard error of an estimated proportion.
# Arguments:
#   p: proportion(s)
#   n: number(s) of datasets (trials)
# Returns: numeric vector of standard errors
mc_se <- function(p, n) sqrt(p * (1 - p) / pmax(n, 1))

# ----------------------------------------------------------------------
# wilson_ci(): Wilson confidence interval for a binomial proportion.
# Arguments:
#   x: number of rejections (successes)
#   n: number(s) of datasets (trials)
#   conf: confidence level
# Returns: numeric vector c(lower, upper)
wilson_ci <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2))) / denom
  c(max(0, centre - half), min(1, centre + half))
}

# ----------------------------------------------------------------------
# simulate_binary_scenario(): Simulate one Bernoulli random-intercept dataset under a head-
#  to-head scenario (null, omitted x2, or serial dependence).
# Arguments:
#   scenario: scenario label
#   m: number of clusters
#   n_i: number of observations per cluster
#   seed: random-number seed
# Returns: data frame with columns id, time, x1, x2, y and scenario
simulate_binary_scenario <- function(scenario, m = M, n_i = N_I, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  id <- rep(seq_len(m), each = n_i)
  time <- rep(seq_len(n_i), times = m)

  x1 <- rnorm(m * n_i)
  x2 <- rnorm(m * n_i)

  b <- rnorm(m, mean = 0, sd = SIGMA_B)

  gamma <- 0
  rho <- 0

  if (scenario == "omitted_x2_gamma05") gamma <- 0.50
  if (scenario == "omitted_x2_gamma10") gamma <- 1.00
  if (scenario == "serial_dependence_rho10") rho <- 1.00

  y <- numeric(m * n_i)

  for (i in seq_len(m)) {
    rows <- which(id == i)
    ylag <- 0

    for (pos in seq_along(rows)) {
      row <- rows[pos]
      eta <- BETA0 + BETA1 * x1[row] + gamma * x2[row] + b[i] + rho * ylag
      p <- plogis(eta)
      y[row] <- rbinom(1, 1, p)
      ylag <- y[row]
    }
  }

  data.frame(
    id = factor(id),
    time = time,
    x1 = x1,
    x2 = x2,
    y = y,
    scenario = scenario
  )
}

# ----------------------------------------------------------------------
# fit_null_binary(): Fit the Bernoulli random-intercept model y ~ x1 + (1 | id) with
#  lme4::glmer.
# Arguments:
#   dat: simulated data frame with columns id, time, x1, x2 and y
# Returns: fitted glmerMod object
fit_null_binary <- function(dat) {
  suppressMessages(
    suppressWarnings(
      glmer(
        y ~ x1 + (1 | id),
        data = dat,
        family = binomial,
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )
    )
  )
}

# ----------------------------------------------------------------------
# rr_binary_glmer_ri(): Sequential posterior-predictive recursive residuals (Eqs. 4-7) for a
#  Bernoulli random-intercept glmer fit, with randomized PIT.
# Arguments:
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   seed: random-number seed
# Returns: data frame (sorted by cluster and order) with added columns rr, u and pred_prob
rr_binary_glmer_ri <- function(fit, data, fixed_formula = ~ x1, seed = 2026) {
  set.seed(seed)
  data <- data[order(data$id, data$time), ]

  beta_hat <- fixef(fit)
  X <- model.matrix(fixed_formula, data = data)
  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]

  vc <- VarCorr(fit)
  sd_b <- as.numeric(attr(vc[["id"]], "stddev"))[1]
  if (!is.finite(sd_b) || sd_b < 1e-8) sd_b <- 1e-8

  gh <- statmod::gauss.quad.prob(n = N_QUAD, dist = "normal")
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights

  y <- as.numeric(data$y)
  ids <- data$id

  rr <- numeric(nrow(data))
  u <- numeric(nrow(data))
  pred_prob <- numeric(nrow(data))

  for (id in unique(ids)) {
    rows <- which(ids == id)
    log_hist <- rep(0, length(b_nodes))

    for (pos in seq_along(rows)) {
      row <- rows[pos]
      logw <- log(base_w) + log_hist
      post_w <- normalize_log_weights(logw)

      eta_nodes <- as.numeric(X[row, , drop = FALSE] %*% beta_hat) + b_nodes
      p_nodes <- plogis(eta_nodes)
      pred_prob[row] <- sum(post_w * p_nodes)

      k <- y[row]
      if (k == 0) {
        G_low <- 0
        G_high <- sum(post_w * (1 - p_nodes))
      } else {
        G_low <- sum(post_w * (1 - p_nodes))
        G_high <- 1
      }

      a <- runif(1)
      u[row] <- G_low + a * (G_high - G_low)
      rr[row] <- safe_qnorm(u[row])
      log_hist <- log_hist + dbinom(k, size = 1, prob = p_nodes, log = TRUE)
    }
  }

  data.frame(data, rr = rr, u = u, pred_prob = pred_prob)
}

# ----------------------------------------------------------------------
# run_rr_calibration(): Refit-calibrated recursive-residual diagnostic for one simulated
#  dataset.
# Arguments:
#   fit: fitted GLMM object
#   dat: simulated data frame with columns id, time, x1, x2 and y
#   rep_id: replicate identifier, used to derive the random-number seeds
# Returns: one-row data frame with rr_pcombo, rr_pmin and rr_best_stat
run_rr_calibration <- function(fit, dat, rep_id) {
  rr_obs <- rr_binary_glmer_ri(fit, dat, seed = 2000000 + rep_id)
  obs_stats <- component_stats(rr_obs, c("time", "pred_prob", "x1", "x2"), "id", "time")

  sim_stats <- vector("list", B_CAL)
  used <- 0

  for (b in seq_len(B_CAL)) {
    one_stats <- tryCatch({
      ysim <- simulate(fit, nsim = 1, seed = 3000000 + 1000 * rep_id + b)[[1]]
      dat_b <- dat
      dat_b$y <- as.numeric(ysim)
      fit_b <- try(fit_null_binary(dat_b), silent = TRUE)
      if (inherits(fit_b, "try-error")) fit_b <- fit
      rr_b <- rr_binary_glmer_ri(fit_b, dat_b, seed = 4000000 + 1000 * rep_id + b)
      component_stats(rr_b, c("time", "pred_prob", "x1", "x2"), "id", "time")
    }, error = function(e) NULL)

    if (!is.null(one_stats)) {
      used <- used + 1
      sim_stats[[used]] <- one_stats
    }
  }

  if (used < max(40, floor(0.5 * B_CAL))) {
    return(data.frame(rr_pcombo = NA_real_, rr_pmin = NA_real_, rr_best_stat = NA_character_))
  }

  sim_stats <- bind_rows(sim_stats[seq_len(used)])
  cal <- calibrate_components(obs_stats, sim_stats)

  data.frame(
    rr_pcombo = cal$p_combo_minP,
    rr_pmin = cal$p_min_raw,
    rr_best_stat = cal$best_stat
  )
}

# ----------------------------------------------------------------------
# run_dharma_tests(): Three non-refit global DHARMa tests used in the manuscript.
#  Temporal autocorrelation is not included in the published comparison.
#  Do not interpret the intentionally missing temporal column as a negative test.
# Arguments:
#   fit: fitted GLMM object
#   dat: simulated data frame with columns id, time, x1, x2 and y
# Returns: one-row data frame of DHARMa p-values and minimum-p summaries
run_dharma_tests <- function(fit, dat) {
  # Original published comparison: three non-refit GLOBAL tests only.
  # A cluster-aware temporal test was NOT successfully computed in the archived
  # study. No temporal p-value enters either combined decision rule.
  out <- list(dharma_uniformity_p = NA_real_, dharma_dispersion_p = NA_real_,
              dharma_outlier_p = NA_real_, dharma_temporal_p = NA_real_,
              dharma_min_raw_p = NA_real_, dharma_min_bonf_p = NA_real_)
  simres <- DHARMa::simulateResiduals(fittedModel = fit, n = DHARMA_NSIM,
                                       refit = FALSE, plot = FALSE)
  out$dharma_uniformity_p <- DHARMa::testUniformity(simres, plot = FALSE)$p.value
  out$dharma_dispersion_p <- DHARMa::testDispersion(simres, plot = FALSE)$p.value
  out$dharma_outlier_p <- DHARMa::testOutliers(simres, plot = FALSE)$p.value
  pvals <- unlist(out[c("dharma_uniformity_p", "dharma_dispersion_p", "dharma_outlier_p")])
  if (length(pvals) != 3L || any(!is.finite(pvals)) || any(pvals < 0 | pvals > 1))
    stop("DHARMa global comparison invalid: all three p-values must be finite in [0,1].")
  out$dharma_min_raw_p <- min(pvals)
  out$dharma_min_bonf_p <- min(1, 3 * min(pvals))
  as.data.frame(out)
}

# ----------------------------------------------------------------------
# run_one_comparison(): Simulate one dataset for a scenario and apply the recursive-residual
#  and DHARMa diagnostics.
# Arguments:
#   scenario: scenario label
#   rep_id: replicate identifier, used to derive the random-number seeds
# Returns: one-row data frame of p-values and detection flags
run_one_comparison <- function(scenario, rep_id) {
  dat <- simulate_binary_scenario(
    scenario = scenario,
    seed = 1000000 + 10000 * match(scenario, SCENARIOS) + rep_id
  )

  fit <- try(fit_null_binary(dat), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(data.frame(
      scenario = scenario,
      rep_id = rep_id,
      fit_failed = TRUE,
      rr_pcombo = NA_real_,
      rr_detect = NA,
      rr_best_stat = NA_character_,
      dharma_uniformity_p = NA_real_,
      dharma_dispersion_p = NA_real_,
      dharma_outlier_p = NA_real_,
      dharma_temporal_p = NA_real_,
      dharma_min_raw_p = NA_real_,
      dharma_min_bonf_p = NA_real_,
      dharma_any_raw_detect = NA,
      dharma_any_bonf_detect = NA
    ))
  }

  rr <- try(run_rr_calibration(fit, dat, rep_id), silent = TRUE)
  if (inherits(rr, "try-error")) {
    rr <- data.frame(rr_pcombo = NA_real_, rr_pmin = NA_real_, rr_best_stat = NA_character_)
  }

  dh <- run_dharma_tests(fit, dat)

  data.frame(
    scenario = scenario,
    rep_id = rep_id,
    fit_failed = FALSE,
    rr_pcombo = rr$rr_pcombo,
    rr_detect = rr$rr_pcombo < ALPHA,
    rr_best_stat = rr$rr_best_stat,
    dharma_uniformity_p = dh$dharma_uniformity_p,
    dharma_dispersion_p = dh$dharma_dispersion_p,
    dharma_outlier_p = dh$dharma_outlier_p,
    dharma_temporal_p = dh$dharma_temporal_p,
    dharma_min_raw_p = dh$dharma_min_raw_p,
    dharma_min_bonf_p = dh$dharma_min_bonf_p,
    dharma_any_raw_detect = dh$dharma_min_raw_p < ALPHA,
    dharma_any_bonf_detect = dh$dharma_min_bonf_p < ALPHA
  )
}

# ----------------------------------------------------------------------
# summarise_detector(): Detection rate, Monte Carlo standard error and Wilson interval of a
#  detector by scenario.
# Arguments:
#   data: data frame used to fit the model
#   detect_var: name of the logical detection column
# Returns: data frame with one row per scenario
summarise_detector <- function(data, detect_var) {
  data |>
    group_by(scenario) |>
    summarise(
      detector = detect_var,
      N_success = sum(!is.na(.data[[detect_var]])),
      detected_n = sum(.data[[detect_var]], na.rm = TRUE),
      detection_rate = detected_n / N_success,
      mc_se = mc_se(detection_rate, N_success),
      ci_low = wilson_ci(detected_n, N_success)[1],
      ci_high = wilson_ci(detected_n, N_success)[2],
      .groups = "drop"
    )
}

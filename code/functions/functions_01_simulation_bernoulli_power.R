# ============================================================
# Function definitions for code/01_simulation_bernoulli_power.R
# Sourced by that script; run it from the root of the supplement.
# ============================================================

# ----------------------------------------------------------------------
# get_output_root(): Return the output folder: RR_GLMM_OUTPUT_ROOT if set, otherwise
#  ./results (created if needed).
# Arguments: none
# Returns: normalized path of the output folder
get_output_root <- function() {
  forced <- Sys.getenv("RR_GLMM_OUTPUT_ROOT", unset = "")
  if (nzchar(forced)) {
    dir.create(forced, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(forced, winslash = "/", mustWork = FALSE))
  }

  default_root <- file.path(getwd(), "results")
  dir.create(default_root, recursive = TRUE, showWarnings = FALSE)
  return(normalizePath(default_root, winslash = "/", mustWork = FALSE))
}

# ----------------------------------------------------------------------
# log_msg(): Write a time-stamped message to the console and to the log file log_file.
# Arguments:
#   ...: message parts, pasted together
# Returns: invisible NULL (called for its side effects)
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
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
    dplyr::group_by(.data[[id_var]]) |>
    dplyr::summarise(C = sum(.data[[rr_var]]) / sqrt(dplyr::n()), .groups = "drop")

  out[["Tmax_C"]] <- max(abs(C$C), na.rm = TRUE)
  out[["Tdisp"]] <- abs(mean(r^2, na.rm = TRUE) - 1)

  rr_lag <- rr_data |>
    dplyr::arrange(.data[[id_var]], .data[[order_var]]) |>
    dplyr::group_by(.data[[id_var]]) |>
    dplyr::mutate(rr_lag = dplyr::lag(.data[[rr_var]])) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(rr_lag))

  out[["Tlag"]] <- if (nrow(rr_lag) > 3) {
    abs(stats::cor(rr_lag[[rr_var]], rr_lag$rr_lag))
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
#   alpha: nominal level for the detection flag
# Returns: data frame with component p-values, p_min, p_combo, leading statistic and
#  detection flag
calibrate_components <- function(obs_stats, sim_stats, alpha = 0.05) {
  stat_names <- intersect(names(obs_stats), names(sim_stats))

  stat_names <- stat_names[
    vapply(obs_stats[stat_names], function(x) is.finite(as.numeric(x[1])), FUN.VALUE = logical(1))
  ]

  stat_names <- stat_names[
    vapply(sim_stats[stat_names], function(x) any(is.finite(x)), FUN.VALUE = logical(1))
  ]

  if (length(stat_names) == 0) {
    return(data.frame(
      statistic = NA_character_,
      observed = NA_real_,
      raw_component_p = NA_real_,
      p_min_raw = NA_real_,
      p_combo_minP = NA_real_,
      best_stat = NA_character_,
      detected_combo = NA
    ))
  }

  raw_p <- vapply(stat_names, function(s) {
    p_from_sim(obs_stats[[s]][1], sim_stats[[s]])
  }, FUN.VALUE = numeric(1))

  p_min_raw <- min(raw_p, na.rm = TRUE)

  sim_p <- matrix(NA_real_, nrow = nrow(sim_stats), ncol = length(stat_names))
  colnames(sim_p) <- stat_names

  for (j in seq_along(stat_names)) {
    s <- stat_names[j]
    vals <- sim_stats[[s]]
    sim_p[, j] <- vapply(vals, function(v) p_from_sim(v, vals), FUN.VALUE = numeric(1))
  }

  sim_min <- apply(sim_p, 1, min, na.rm = TRUE)
  p_combo <- (1 + sum(sim_min <= p_min_raw, na.rm = TRUE)) / (length(sim_min) + 1)

  data.frame(
    statistic = stat_names,
    observed = as.numeric(obs_stats[1, stat_names]),
    raw_component_p = as.numeric(raw_p),
    p_min_raw = p_min_raw,
    p_combo_minP = p_combo,
    best_stat = names(which.min(raw_p)),
    detected_combo = p_combo < alpha
  )
}

# ----------------------------------------------------------------------
# mode_string(): Most frequent non-missing value of a character vector.
# Arguments:
#   x: number of rejections (successes)
# Returns: the most frequent value, or NA
mode_string <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

# ----------------------------------------------------------------------
# mc_se(): Monte Carlo standard error of an estimated proportion.
# Arguments:
#   p: proportion(s)
#   n: number(s) of datasets (trials)
# Returns: numeric vector of standard errors
mc_se <- function(p, n) {
  ifelse(is.finite(p) & is.finite(n) & n > 0, sqrt(p * (1 - p) / n), NA_real_)
}

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
# format_pct_ci(): Format a proportion and its interval as a percentage string.
# Arguments:
#   p: proportion(s)
#   lo: lower interval limit (proportion)
#   hi: upper interval limit (proportion)
# Returns: character string, e.g. "6.2% (5.6, 7.0)"
format_pct_ci <- function(p, lo, hi) {
  if (!is.finite(p) || !is.finite(lo) || !is.finite(hi)) return(NA_character_)
  sprintf("%.1f%% (%.1f, %.1f)", 100 * p, 100 * lo, 100 * hi)
}

# ----------------------------------------------------------------------
# simulate_bernoulli_ri(): Simulate a Bernoulli random-intercept dataset with covariates x1
#  and x2.
# Arguments:
#   m: number of clusters
#   n_i: number of observations per cluster
#   gamma: coefficient of the omitted covariate x2 in the data-generating model
#   beta0: intercept of the data-generating model
#   beta1: coefficient of x1 in the data-generating model
#   sigma_b: standard deviation of the random intercept
#   seed: random-number seed
# Returns: data frame with columns id, time, x1, x2, y, true_p, gamma, m and n_i
simulate_bernoulli_ri <- function(m, n_i, gamma,
                                  beta0 = BETA0,
                                  beta1 = BETA1,
                                  sigma_b = SIGMA_B,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  id <- rep(seq_len(m), each = n_i)
  time <- rep(seq_len(n_i), times = m)

  x1 <- rnorm(m * n_i)
  x2 <- rnorm(m * n_i)

  b <- rnorm(m, mean = 0, sd = sigma_b)

  eta <- beta0 + beta1 * x1 + gamma * x2 + b[id]
  p <- plogis(eta)
  y <- rbinom(m * n_i, size = 1, prob = p)

  data.frame(
    id = factor(id),
    time = time,
    x1 = x1,
    x2 = x2,
    y = y,
    true_p = p,
    gamma = gamma,
    m = m,
    n_i = n_i
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
      lme4::glmer(
        y ~ x1 + (1 | id),
        data = dat,
        family = stats::binomial,
        control = lme4::glmerControl(
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
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   n_quad: number of Gauss-Hermite quadrature nodes
#   seed: random-number seed
# Returns: data frame (sorted by cluster and order) with added columns rr, u and pred_prob
rr_binary_glmer_ri <- function(fit, data,
                               fixed_formula = ~ x1,
                               id_var = "id",
                               order_var = "time",
                               y_var = "y",
                               n_quad = N_QUAD,
                               seed = 2026) {
  set.seed(seed)

  data <- data[order(data[[id_var]], data[[order_var]]), ]

  beta_hat <- lme4::fixef(fit)
  X <- stats::model.matrix(fixed_formula, data = data)

  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]

  vc <- lme4::VarCorr(fit)
  sd_b <- as.numeric(attr(vc[[id_var]], "stddev"))[1]
  if (!is.finite(sd_b) || sd_b < 1e-8) sd_b <- 1e-8

  gh <- statmod::gauss.quad.prob(n = n_quad, dist = "normal")
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights

  y <- as.numeric(data[[y_var]])
  ids <- data[[id_var]]

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

      a <- stats::runif(1)
      u[row] <- G_low + a * (G_high - G_low)
      rr[row] <- safe_qnorm(u[row])

      log_hist <- log_hist + stats::dbinom(k, size = 1, prob = p_nodes, log = TRUE)
    }
  }

  data.frame(data, rr = rr, u = u, pred_prob = pred_prob)
}

# ----------------------------------------------------------------------
# run_one_replicate(): Simulate one dataset, fit the null model and compute the refit-
#  calibrated omnibus test.
# Arguments:
#   m: number of clusters
#   n_i: number of observations per cluster
#   gamma: coefficient of the omitted covariate x2 in the data-generating model
#   rep_id: replicate identifier, used to derive the random-number seeds
#   row_type: label stored in the output ("null_level" or "power")
# Returns: one-row data frame with design settings, calibration counts, p_combo, p_min and
#  leading statistic
run_one_replicate <- function(m, n_i, gamma, rep_id, row_type) {
  dat <- simulate_bernoulli_ri(
    m = m,
    n_i = n_i,
    gamma = gamma,
    seed = 1000000 + 10000 * m + 100 * n_i + rep_id + round(1000 * gamma)
  )

  fit <- try(fit_null_binary(dat), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(data.frame(
      row_type = row_type,
      m = m, n_i = n_i, gamma = gamma, rep_id = rep_id,
      fit_failed = TRUE,
      calibration_successful = NA_integer_,
      refit_failures = NA_integer_,
      rr_failures = NA_integer_,
      singular_refits = NA_integer_,
      p_combo_minP = NA_real_,
      p_min_raw = NA_real_,
      best_stat = NA_character_,
      detected_combo = NA
    ))
  }

  index_vars <- c("time", "pred_prob", "x1", "x2")

  rr_obs <- try(rr_binary_glmer_ri(fit, dat, seed = 2000000 + rep_id), silent = TRUE)

  if (inherits(rr_obs, "try-error")) {
    return(data.frame(
      row_type = row_type,
      m = m, n_i = n_i, gamma = gamma, rep_id = rep_id,
      fit_failed = FALSE,
      calibration_successful = 0L,
      refit_failures = NA_integer_,
      rr_failures = NA_integer_,
      singular_refits = NA_integer_,
      p_combo_minP = NA_real_,
      p_min_raw = NA_real_,
      best_stat = NA_character_,
      detected_combo = NA
    ))
  }

  obs_stats <- component_stats(rr_obs, index_vars, "id", "time")
  sim_stats <- vector("list", B_CAL)

  used <- 0
  refit_fail <- 0
  rr_fail <- 0
  singular_count <- 0

  for (b in seq_len(B_CAL)) {
    one_stats <- tryCatch({
      ysim <- stats::simulate(fit, nsim = 1, seed = 3000000 + 1000 * rep_id + b)[[1]]

      dat_b <- dat
      dat_b$y <- as.numeric(ysim)

      fit_b <- try(fit_null_binary(dat_b), silent = TRUE)
      if (inherits(fit_b, "try-error")) {
        refit_fail <<- refit_fail + 1
        fit_b <- fit
      }

      if (lme4::isSingular(fit_b)) singular_count <<- singular_count + 1

      rr_b <- rr_binary_glmer_ri(
        fit_b,
        dat_b,
        seed = 4000000 + 1000 * rep_id + b
      )

      component_stats(rr_b, index_vars, "id", "time")
    }, error = function(e) {
      rr_fail <<- rr_fail + 1
      NULL
    })

    if (!is.null(one_stats)) {
      used <- used + 1
      sim_stats[[used]] <- one_stats
    }
  }

  if (used < max(40, floor(0.50 * B_CAL))) {
    return(data.frame(
      row_type = row_type,
      m = m, n_i = n_i, gamma = gamma, rep_id = rep_id,
      fit_failed = FALSE,
      calibration_successful = used,
      refit_failures = refit_fail,
      rr_failures = rr_fail,
      singular_refits = singular_count,
      p_combo_minP = NA_real_,
      p_min_raw = NA_real_,
      best_stat = NA_character_,
      detected_combo = NA
    ))
  }

  sim_stats <- dplyr::bind_rows(sim_stats[seq_len(used)])
  rr_cal <- calibrate_components(obs_stats, sim_stats, ALPHA)

  data.frame(
    row_type = row_type,
    m = m,
    n_i = n_i,
    gamma = gamma,
    rep_id = rep_id,
    fit_failed = FALSE,
    calibration_successful = used,
    refit_failures = refit_fail,
    rr_failures = rr_fail,
    singular_refits = singular_count,
    p_combo_minP = rr_cal$p_combo_minP[1],
    p_min_raw = rr_cal$p_min_raw[1],
    best_stat = rr_cal$best_stat[1],
    detected_combo = rr_cal$detected_combo[1]
  )
}

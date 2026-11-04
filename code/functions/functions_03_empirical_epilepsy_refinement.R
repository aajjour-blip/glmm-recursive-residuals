# ============================================================
# Function definitions for code/03_empirical_epilepsy_refinement.R
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
  normalizePath(default_root, winslash = "/", mustWork = FALSE)
}

# ----------------------------------------------------------------------
# normalize_log_weights(): Convert log weights to normalized weights, stabilized by
#  subtracting the maximum.
# Arguments:
#   logw: numeric vector of log weights at the quadrature nodes
# Returns: numeric vector of weights summing to one
normalize_log_weights <- function(logw) { m <- max(logw); w <- exp(logw - m); w / sum(w) }

# ----------------------------------------------------------------------
# safe_qnorm(): Standard normal quantile function with probabilities truncated to [eps, 1 -
#  eps].
# Arguments:
#   u: numeric vector of probabilities
# Returns: numeric vector of normal quantiles
safe_qnorm <- function(u) { eps <- .Machine$double.eps; qnorm(pmin(pmax(u, eps), 1 - eps)) }

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
  data.frame(q = rr_data[[q_var]][ord], W = cumsum(rr_data[[rr_var]][ord]) / sqrt(nrow(rr_data)))
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
    if (all(is.na(rr_data[[v]])) || length(unique(rr_data[[v]])) < 2) next
    proc <- rr_process(rr_data, v, rr_var)
    nm <- make.names(v)
    out[[paste0("Tinf_", nm)]] <- max(abs(proc$W), na.rm = TRUE)
    out[[paste0("T2_", nm)]] <- mean(proc$W^2, na.rm = TRUE)
  }
  r <- rr_data[[rr_var]]
  out[["Tmax_r"]] <- max(abs(r), na.rm = TRUE)
  C <- rr_data |> group_by(.data[[id_var]]) |> summarise(C = sum(.data[[rr_var]]) / sqrt(n()),
                                                         .groups = "drop")
  out[["Tmax_C"]] <- max(abs(C$C), na.rm = TRUE)
  out[["Tdisp"]] <- abs(mean(r^2, na.rm = TRUE) - 1)
  rr_lag <- rr_data |> arrange(.data[[id_var]], .data[[order_var]]) |>
    group_by(.data[[id_var]]) |> mutate(rr_lag = lag(.data[[rr_var]])) |>
    ungroup() |> filter(!is.na(rr_lag))
  out[["Tlag"]] <- if (nrow(rr_lag) > 3) abs(cor(rr_lag[[rr_var]], rr_lag$rr_lag)) else NA_real_
  as.data.frame(out)
}

# ----------------------------------------------------------------------
# p_from_sim(): Refit-calibrated upper-tail p-value of an observed statistic, (1 + #{sim >=
#  obs}) / (B + 1).
# Arguments:
#   obs: observed value of a statistic
#   sim: numeric vector of reference values of the same statistic
# Returns: p-value, or NA if no finite reference values
p_from_sim <- function(obs, sim) (1 + sum(sim >= obs,
                                          na.rm = TRUE)) / (length(sim[is.finite(sim)]) + 1)

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
  stat_names <- stat_names[sapply(obs_stats[stat_names], function(x) is.finite(as.numeric(x[1])))]
  stat_names <- stat_names[sapply(sim_stats[stat_names], function(x) any(is.finite(x)))]
  raw_p <- sapply(stat_names, function(s) p_from_sim(obs_stats[[s]][1], sim_stats[[s]]))
  p_min_raw <- min(raw_p, na.rm = TRUE)
  sim_p <- matrix(NA_real_, nrow = nrow(sim_stats), ncol = length(stat_names))
  colnames(sim_p) <- stat_names
  for (j in seq_along(stat_names)) {
    vals <- sim_stats[[stat_names[j]]]
    sim_p[, j] <- sapply(vals, function(v) p_from_sim(v, vals))
  }
  sim_min <- apply(sim_p, 1, min, na.rm = TRUE)
  p_combo <- (1 + sum(sim_min <= p_min_raw, na.rm = TRUE)) / (length(sim_min) + 1)
  data.frame(statistic = stat_names,
             observed = as.numeric(obs_stats[1, stat_names]),
             raw_component_p = as.numeric(raw_p),
             p_min_raw = p_min_raw,
             p_combo_minP = p_combo,
             best_stat = names(which.min(raw_p)),
             detected_combo = p_combo < alpha)
}

# ----------------------------------------------------------------------
# plot_rr_process(): Plot the cumulative residual process for one index and save it as a PNG
#  file.
# Arguments:
#   rr_data: data frame of recursive residuals and index variables (one row per observation)
#   q_var: name of the column used as the index of the cumulative process
#   out_file: path of the PNG file to write
#   title: plot title
# Returns: invisible (called for its side effect)
plot_rr_process <- function(rr_data, q_var, out_file, title) {
  proc <- rr_process(rr_data, q_var)
  p <- ggplot(proc, aes(x = q, y = W)) + geom_line() + geom_point(size = 1) +
    geom_hline(yintercept = 0, linetype = 2) + theme_bw() +
    labs(x = q_var, y = "Recursive residual process", title = title)
  ggsave(out_file, p, width = 7, height = 4.5, dpi = 300)
}

# ----------------------------------------------------------------------
# safe_test(): Evaluate a test expression and return NA if it fails.
# Arguments:
#   expr: expression to evaluate
# Returns: value of expr, or NA
safe_test <- function(expr) tryCatch(expr, error = function(e) NA_real_)

# ----------------------------------------------------------------------
# get_cond_re_sd(): Random-intercept standard deviation of a glmmTMB fit, bounded below by
#  1e-8.
# Arguments:
#   fit: fitted GLMM object
#   id_var: name of the cluster identifier column
# Returns: numeric scalar
get_cond_re_sd <- function(fit, id_var) {
  vc <- VarCorr(fit)$cond[[id_var]]
  sd_b <- as.numeric(attr(vc, "stddev"))[1]
  if (!is.finite(sd_b) || sd_b < 1e-8) sd_b <- 1e-8
  sd_b
}

# ----------------------------------------------------------------------
# get_nb2_theta(): NB2 size parameter of a glmmTMB fit (sigma(fit)).
# Arguments:
#   fit: fitted GLMM object
# Returns: numeric scalar
get_nb2_theta <- function(fit) {
  theta <- tryCatch(as.numeric(sigma(fit)), error = function(e) NA_real_)
  if (!is.finite(theta) || theta <= 0) stop("Could not extract valid NB2 theta via sigma(fit).")
  theta
}

# ----------------------------------------------------------------------
# count_pmf(): Probability mass function of a (zero-inflated) Poisson or NB2 count
#  distribution.
# Arguments:
#   k: observed count(s)
#   lambda: conditional mean(s) of the count distribution
#   family_type: count family, "poisson" or "nbinom2"
#   theta: NB2 size (dispersion) parameter; ignored for Poisson
#   zi: zero-inflation probability (0 for no zero inflation)
# Returns: numeric vector of probabilities
count_pmf <- function(k, lambda, family_type, theta = NULL, zi = 0) {
  base <- if (family_type == "poisson") dpois(k, lambda = lambda) else dnbinom(k, mu = lambda,
      size = theta)
  if (length(zi) == 1) zi <- rep(zi, length(lambda))
  if (k == 0) zi + (1 - zi) * base else (1 - zi) * base
}

# ----------------------------------------------------------------------
# count_cdf(): Distribution function of a (zero-inflated) Poisson or NB2 count distribution.
# Arguments:
#   k: observed count(s)
#   lambda: conditional mean(s) of the count distribution
#   family_type: count family, "poisson" or "nbinom2"
#   theta: NB2 size (dispersion) parameter; ignored for Poisson
#   zi: zero-inflation probability (0 for no zero inflation)
# Returns: numeric vector of probabilities
count_cdf <- function(k, lambda, family_type, theta = NULL, zi = 0) {
  if (k < 0) return(rep(0, length(lambda)))
  base <- if (family_type == "poisson") ppois(k, lambda = lambda) else pnbinom(k, mu = lambda,
      size = theta)
  if (length(zi) == 1) zi <- rep(zi, length(lambda))
  zi + (1 - zi) * base
}

# ----------------------------------------------------------------------
# rr_count_glmmTMB_ri(): Sequential posterior-predictive recursive residuals for a Poisson,
#  NB2 or zero-inflated random-intercept glmmTMB fit, with randomized
#  PIT.
# Arguments:
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   family_type: count family, "poisson" or "nbinom2"
#   zero_inflated: logical; TRUE for a zero-inflated count model
#   zi_formula: one-sided formula defining the zero-inflation design matrix
#   offset_fun: function of the data returning offset values, or NULL
#   n_quad: number of Gauss-Hermite quadrature nodes
#   seed: random-number seed
# Returns: data frame (sorted by cluster and order) with added columns rr, u, pred_mean and
#  zi_prob
rr_count_glmmTMB_ri <- function(fit, data, fixed_formula, id_var, order_var, y_var,
                                family_type = c("poisson", "nbinom2"), zero_inflated = FALSE,
                                zi_formula = NULL, offset_fun = NULL,
                                n_quad = 40, seed = 2026) {
  family_type <- match.arg(family_type)
  set.seed(seed)
  data <- data[order(data[[id_var]], data[[order_var]]), ]

  beta_hat <- fixef(fit)$cond
  X <- model.matrix(fixed_formula, data)
  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]
  missing_beta <- setdiff(names(fixef(fit)$cond), colnames(X))
  if (length(missing_beta) > 0) stop("Missing fixed-effect columns in conditional X: ",
                                     paste(missing_beta, collapse = ", "))

  if (zero_inflated) {
    beta_zi <- fixef(fit)$zi
    if (is.null(zi_formula)) stop("zi_formula is required for zero-inflated RR.")
    X_zi <- model.matrix(zi_formula, data)
    common_zi <- intersect(colnames(X_zi), names(beta_zi))
    X_zi <- X_zi[, common_zi, drop = FALSE]
    beta_zi <- beta_zi[common_zi]
    missing_zi <- setdiff(names(fixef(fit)$zi), colnames(X_zi))
    if (length(missing_zi) > 0) stop("Missing fixed-effect columns in ZI X: ",
                                     paste(missing_zi, collapse = ", "))
    zi_prob_all <- as.numeric(plogis(X_zi %*% beta_zi))
  } else {
    zi_prob_all <- rep(0, nrow(data))
  }

  offset_values <- if (is.null(offset_fun)) rep(0, nrow(data)) else offset_fun(data)
  theta <- if (family_type == "nbinom2") get_nb2_theta(fit) else NULL
  sd_b <- get_cond_re_sd(fit, id_var)

  gh <- statmod::gauss.quad.prob(n = n_quad, dist = "normal")
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights
  y <- data[[y_var]]
  ids <- data[[id_var]]

  rr <- numeric(nrow(data)); u <- numeric(nrow(data)); pred_mean <- numeric(nrow(data))

  for (id in unique(ids)) {
    rows <- which(ids == id)
    log_hist <- rep(0, length(b_nodes))
    for (pos in seq_along(rows)) {
      row <- rows[pos]
      post_w <- normalize_log_weights(log(base_w) + log_hist)
      eta_nodes <- as.numeric(X[row, , drop = FALSE] %*% beta_hat) + b_nodes + offset_values[row]
      lambda_nodes <- exp(eta_nodes)
      zi_nodes <- rep(zi_prob_all[row], length(lambda_nodes))
      k <- y[row]
      G_low <- sum(post_w * count_cdf(k - 1, lambda_nodes, family_type, theta, zi_nodes))
      G_high <- sum(post_w * count_cdf(k, lambda_nodes, family_type, theta, zi_nodes))
      a <- runif(1)
      u[row] <- G_low + a * (G_high - G_low)
      rr[row] <- safe_qnorm(u[row])
      pred_mean[row] <- sum(post_w * (1 - zi_nodes) * lambda_nodes)
      pmf_nodes <- pmax(count_pmf(k, lambda_nodes, family_type, theta, zi_nodes),
                        .Machine$double.xmin)
      log_hist <- log_hist + log(pmf_nodes)
    }
  }
  data.frame(data, rr = rr, u = u, pred_mean = pred_mean, zi_prob = zi_prob_all)
}

# ----------------------------------------------------------------------
# run_count_rr_model(): Observed recursive residuals, refit calibration and output files for
#  one count model.
# Arguments:
#   example: example label used in output file names
#   model_name: model label used in output file names
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   refit_fun: function that refits the model to a data frame with a simulated response
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   family_type: count family, "poisson" or "nbinom2"
#   zero_inflated: logical; TRUE for a zero-inflated count model
#   zi_formula: one-sided formula defining the zero-inflation design matrix
#   offset_fun: function of the data returning offset values, or NULL
#   index_vars: character vector of index columns for the cumulative processes
# Returns: one-row data frame summarizing the run (p_combo, p_min, leading statistic,
#  failures)
run_count_rr_model <- function(example, model_name, fit, data, fixed_formula, refit_fun,
                               id_var, order_var, y_var, family_type, zero_inflated = FALSE,
                               zi_formula = NULL, offset_fun = NULL, index_vars) {
  cat("\n============================================================\n")
  cat("RR-guided refinement:", example, model_name, "\n")
  cat("============================================================\n")

  rr_obs <- rr_count_glmmTMB_ri(fit, data, fixed_formula, id_var, order_var, y_var,
                                family_type, zero_inflated, zi_formula, offset_fun,
                                N_QUAD, 2026)
  obs_stats <- component_stats(rr_obs, index_vars, id_var, order_var)
  sim_stats <- vector("list", B_ENV)
  refit_fail <- 0; rr_fail <- 0; used <- 0

  for (b in seq_len(B_ENV)) {
    if (b %% 50 == 0) cat(example, model_name, "RR simulation", b, "of", B_ENV, "\n")
    one_stats <- tryCatch({
      ysim <- simulate(fit, nsim = 1, seed = 90000 + b)[[1]]
      dat_b <- data
      dat_b[[y_var]] <- as.numeric(ysim)
      fit_b <- try(suppressWarnings(refit_fun(dat_b)), silent = TRUE)
      if (inherits(fit_b, "try-error")) { refit_fail <<- refit_fail + 1; fit_b <- fit }
      rr_b <- rr_count_glmmTMB_ri(fit_b, dat_b, fixed_formula, id_var, order_var, y_var,
                                  family_type, zero_inflated, zi_formula, offset_fun,
                                  N_QUAD, 100000 + b)
      component_stats(rr_b, index_vars, id_var, order_var)
    }, error = function(e) { rr_fail <<- rr_fail + 1; NULL })
    if (!is.null(one_stats)) { used <- used + 1; sim_stats[[used]] <- one_stats }
  }

  sim_stats <- bind_rows(sim_stats[seq_len(used)])
  if (nrow(sim_stats) < 50) stop(example, " ", model_name,
                                 ": too few successful calibration replicates: ", nrow(sim_stats))
  rr_cal <- calibrate_components(obs_stats, sim_stats, ALPHA)
  prefix <- paste0(example, "_", model_name)
  write_csv(rr_cal, file.path(OUT, "data", paste0(prefix, "_RR_B1000_component_calibration.csv")))
  write_csv(rr_obs, file.path(OUT, "data", paste0(prefix, "_RR_observed_residuals.csv")))
  write_csv(sim_stats, file.path(OUT, "data", paste0(prefix, "_RR_B1000_simulated_stats.csv")))

  for (v in index_vars[seq_len(min(5, length(index_vars)))]) {
    if (v %in% names(rr_obs)) plot_rr_process(rr_obs, v, file.path(OUT, "figures",
        paste0(prefix, "_RR_process_", make.names(v), ".png")), paste(prefix, "RR process by", v))
  }

  summary <- data.frame(example = example, model_name = model_name, family_type = family_type,
                        zero_inflated = zero_inflated, B_ENV_requested = B_ENV,
                        B_ENV_successful = nrow(sim_stats), N_QUAD = N_QUAD,
                        refit_failures = refit_fail, rr_failures_skipped = rr_fail,
                        p_combo_minP = rr_cal$p_combo_minP[1], p_min_raw = rr_cal$p_min_raw[1],
                        best_stat = rr_cal$best_stat[1], detected_combo = rr_cal$detected_combo[1])
  write_csv(summary, file.path(OUT, "data", paste0(prefix, "_RR_run_summary.csv")))
  summary
}

# ----------------------------------------------------------------------
# dharma_for_model(): DHARMa residual tests for one model; also saves the DHARMa plot.
# Arguments:
#   fit: fitted GLMM object
#   example: example label used in output file names
#   model_name: model label used in output file names
# Returns: one-row data frame of DHARMa p-values
dharma_for_model <- function(fit, example, model_name) {
  cat("DHARMa for", example, model_name, "\n")
  sim <- DHARMa::simulateResiduals(fit, n = B_ENV, seed = 2026)
  png(file.path(OUT, "figures", paste0(example, "_", model_name, "_DHARMa.png")), 1200, 900)
  plot(sim); dev.off()
  data.frame(example = example, model_name = model_name,
             DHARMa_uniformity_p = safe_test(DHARMa::testUniformity(sim)$p.value),
             DHARMa_dispersion_p = safe_test(DHARMa::testDispersion(sim)$p.value),
             DHARMa_zeroInflation_p = safe_test(DHARMa::testZeroInflation(sim)$p.value),
             DHARMa_outlier_p = safe_test(DHARMa::testOutliers(sim, type = "bootstrap")$p.value))
}

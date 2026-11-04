# ============================================================
# Function definitions for code/05_empirical_additional_examples.R
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
  (1 + sum(sim >= obs, na.rm = TRUE)) / (length(sim[is.finite(sim)]) + 1)
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
  stat_names <- stat_names[sapply(obs_stats[stat_names], function(x) is.finite(as.numeric(x[1])))]
  stat_names <- stat_names[sapply(sim_stats[stat_names], function(x) any(is.finite(x)))]

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
  p <- ggplot2::ggplot(proc, ggplot2::aes(x = q, y = W)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = q_var, y = "Recursive residual process", title = title)

  ggplot2::ggsave(out_file, p, width = 7, height = 4.5, dpi = 300)
}

# ----------------------------------------------------------------------
# safe_test(): Evaluate a test expression and return NA if it fails.
# Arguments:
#   expr: expression to evaluate
# Returns: value of expr, or NA
safe_test <- function(expr) {
  tryCatch(expr, error = function(e) NA_real_)
}

# ----------------------------------------------------------------------
# write_model_AIC(): Compute AIC and BIC for a list of models and write them to a CSV file.
# Arguments:
#   models: named list of fitted models
#   file: path of the output CSV file
# Returns: data frame of information criteria
write_model_AIC <- function(models, file) {
  # current version:
  # Do not pass a named list directly to AIC(), because some R versions
  # dispatch to AIC.default() and complain that argument "object" is missing.
  # Instead pass an unnamed list, then restore the model labels explicitly.
  aic_raw <- do.call(stats::AIC, unname(models)) |>
    as.data.frame()

  aic <- data.frame(
    model = names(models),
    df = aic_raw$df,
    AIC = aic_raw$AIC,
    row.names = NULL
  )

  readr::write_csv(aic, file)
  aic
}

# ----------------------------------------------------------------------
# dharma_summary(): DHARMa tests and plots for a named list of models; writes a CSV file.
# Arguments:
#   named_models: named list of fitted models
#   n: number(s) of datasets (trials)
#   prefix: prefix used in output file names
# Returns: data frame of DHARMa p-values, one row per model
dharma_summary <- function(named_models, n = B_ENV, prefix) {
  out <- list()

  for (nm in names(named_models)) {
    cat("DHARMa:", prefix, nm, "\n")
    fit <- named_models[[nm]]

    sim <- DHARMa::simulateResiduals(fit, n = n, seed = 3000 + match(nm, names(named_models)))

    png(file.path(OUT, "figures", paste0(prefix, "_DHARMa_", nm, ".png")), 1200, 900)
    plot(sim)
    dev.off()

    out[[nm]] <- data.frame(
      model = nm,
      DHARMa_uniformity_p = safe_test(DHARMa::testUniformity(sim)$p.value),
      DHARMa_dispersion_p = safe_test(DHARMa::testDispersion(sim)$p.value),
      DHARMa_zeroInflation_p = safe_test(DHARMa::testZeroInflation(sim)$p.value),
      DHARMa_outlier_p = safe_test(DHARMa::testOutliers(sim, type = "bootstrap")$p.value)
    )
  }

  res <- dplyr::bind_rows(out)
  readr::write_csv(res, file.path(OUT, "data", paste0(prefix, "_DHARMa_tests_B1000.csv")))
  res
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
rr_binary_glmer_ri <- function(fit, data, fixed_formula,
                               id_var, order_var, y_var,
                               n_quad = 40, seed = 2026) {
  set.seed(seed)

  data <- data[order(data[[id_var]], data[[order_var]]), ]

  beta_hat <- lme4::fixef(fit)
  X <- stats::model.matrix(fixed_formula, data = data)

  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]

  missing_beta <- setdiff(names(lme4::fixef(fit)), colnames(X))
  if (length(missing_beta) > 0) {
    stop("Missing fixed-effect columns in binary X: ", paste(missing_beta, collapse = ", "))
  }

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
      p_nodes <- stats::plogis(eta_nodes)

      pi_hat <- sum(post_w * p_nodes)
      pred_prob[row] <- pi_hat

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
# run_binary_lme4_rr(): Fit a Bernoulli random-intercept model, run DHARMa, and compute the
#  refit-calibrated recursive-residual diagnostic.
# Arguments:
#   dataset_name: dataset label used in output file names
#   data: data frame used to fit the model
#   formula: model formula passed to lme4::glmer
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   index_vars: character vector of index columns for the cumulative processes
# Returns: one-row data frame summarizing the run
run_binary_lme4_rr <- function(dataset_name, data, formula, fixed_formula,
                               id_var, order_var, y_var, index_vars, diagnostic_queue) {
  cat("\n============================================================\n")
  cat("Binary RR example:", dataset_name, "\n")
  cat("============================================================\n")

  fit <- lme4::glmer(formula, data = data, family = binomial,
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 2e5)))

  # Register comparator work, but do not execute it until every RR calibration
  # in the parent workflow has completed. This avoids DHARMa-version-dependent
  # changes in the random-number stream consumed by later examples.
  stopifnot(is.environment(diagnostic_queue),
            is.list(diagnostic_queue$tasks))
  diagnostic_queue$tasks[[length(diagnostic_queue$tasks) + 1L]] <- local({
    fit_dh <- fit
    dataset_dh <- dataset_name
    formula_dh <- formula
    function() {
      sim <- DHARMa::simulateResiduals(fit_dh, n = B_ENV, seed = 2026)
      png(file.path(OUT, "figures", paste0(dataset_dh, "_DHARMa.png")),
          1200, 900)
      plot(sim)
      dev.off()
      dh <- data.frame(
        dataset = dataset_dh,
        model = deparse(formula_dh),
        DHARMa_uniformity_p = safe_test(DHARMa::testUniformity(sim)$p.value),
        DHARMa_dispersion_p = safe_test(DHARMa::testDispersion(sim)$p.value),
        DHARMa_outlier_p = safe_test(DHARMa::testOutliers(sim, type = "bootstrap")$p.value),
        AIC = AIC(fit_dh),
        BIC = BIC(fit_dh),
        isSingular = lme4::isSingular(fit_dh)
      )
      readr::write_csv(
        dh, file.path(OUT, "data", paste0(dataset_dh, "_DHARMa_AIC_B1000.csv"))
      )
    }
  })

  rr_obs <- rr_binary_glmer_ri(fit, data, fixed_formula,
                               id_var, order_var, y_var, N_QUAD, seed = 2026)
  obs_stats <- component_stats(rr_obs, index_vars, id_var, order_var)

  sim_stats <- vector("list", B_ENV)
  refit_fail <- 0
  rr_fail <- 0
  singular_count <- 0
  used <- 0

  for (b in seq_len(B_ENV)) {
    if (b %% 50 == 0) cat(dataset_name, "RR simulation", b, "of", B_ENV, "\n")

    one_stats <- tryCatch({
      ysim <- simulate(fit, nsim = 1, seed = 10000 + b)[[1]]
      dat_b <- data
      dat_b[[y_var]] <- as.numeric(ysim)

      fit_b <- try(
        suppressWarnings(
          lme4::glmer(formula, data = dat_b, family = binomial,
                      control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
        ),
        silent = TRUE
      )

      if (inherits(fit_b, "try-error")) {
        refit_fail <<- refit_fail + 1
        fit_b <- fit
      }

      if (lme4::isSingular(fit_b)) singular_count <<- singular_count + 1

      rr_b <- rr_binary_glmer_ri(fit_b, dat_b, fixed_formula,
                                 id_var, order_var, y_var,
                                 N_QUAD, seed = 20000 + b)

      component_stats(rr_b, index_vars, id_var, order_var)
    }, error = function(e) {
      rr_fail <<- rr_fail + 1
      NULL
    })

    if (!is.null(one_stats)) {
      used <- used + 1
      sim_stats[[used]] <- one_stats
    }
  }

  sim_stats <- sim_stats[seq_len(used)]
  sim_stats <- dplyr::bind_rows(sim_stats)

  if (nrow(sim_stats) < 50) {
    stop(dataset_name, ": too few successful simulated calibration replicates: ", nrow(sim_stats))
  }

  rr_cal <- calibrate_components(obs_stats, sim_stats, ALPHA)

  readr::write_csv(rr_cal, file.path(OUT, "data", paste0(dataset_name,
                                                         "_RR_B1000_component_calibration.csv")))
  readr::write_csv(rr_obs, file.path(OUT, "data", paste0(dataset_name,
                                                         "_RR_observed_residuals.csv")))
  readr::write_csv(sim_stats, file.path(OUT, "data", paste0(dataset_name,
                                                            "_RR_B1000_simulated_stats.csv")))

  for (v in index_vars[seq_len(min(3, length(index_vars)))]) {
    if (v %in% names(rr_obs)) {
      plot_rr_process(rr_obs, v,
                      file.path(OUT, "figures", paste0(dataset_name, "_RR_process_",
                                                       make.names(v), ".png")),
                      paste(dataset_name, "RR process by", v))
    }
  }

  summary <- data.frame(
    dataset = dataset_name,
    family = "binomial",
    B_ENV_requested = B_ENV,
    B_ENV_successful = nrow(sim_stats),
    N_QUAD = N_QUAD,
    refit_failures = refit_fail,
    rr_failures_skipped = rr_fail,
    singular_refits = singular_count,
    p_combo_minP = rr_cal$p_combo_minP[1],
    p_min_raw = rr_cal$p_min_raw[1],
    best_stat = rr_cal$best_stat[1],
    detected_combo = rr_cal$detected_combo[1]
  )

  readr::write_csv(summary, file.path(OUT, "data", paste0(dataset_name, "_RR_run_summary.csv")))
  summary
}

# ----------------------------------------------------------------------
# get_cond_re_sd(): Random-intercept standard deviation of a glmmTMB fit, bounded below by
#  1e-8.
# Arguments:
#   fit: fitted GLMM object
#   id_var: name of the cluster identifier column
# Returns: numeric scalar
get_cond_re_sd <- function(fit, id_var) {
  vc <- glmmTMB::VarCorr(fit)$cond[[id_var]]
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
  if (!is.finite(theta) || theta <= 0) {
    stop("Could not extract valid NB2 theta via sigma(fit).")
  }
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
  if (family_type == "poisson") {
    base <- stats::dpois(k, lambda = lambda)
  } else if (family_type == "nbinom2") {
    base <- stats::dnbinom(k, mu = lambda, size = theta)
  } else {
    stop("Unsupported family_type: ", family_type)
  }

  if (length(zi) == 1) zi <- rep(zi, length(lambda))

  if (k == 0) {
    zi + (1 - zi) * base
  } else {
    (1 - zi) * base
  }
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

  if (family_type == "poisson") {
    base <- stats::ppois(k, lambda = lambda)
  } else if (family_type == "nbinom2") {
    base <- stats::pnbinom(k, mu = lambda, size = theta)
  } else {
    stop("Unsupported family_type: ", family_type)
  }

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
#   offset_values: numeric vector of offsets in the row order of data, or NULL
#   n_quad: number of Gauss-Hermite quadrature nodes
#   seed: random-number seed
# Returns: data frame (sorted by cluster and order) with added columns rr, u, pred_mean and
#  zi_prob
rr_count_glmmTMB_ri <- function(fit, data, fixed_formula,
                                id_var, order_var, y_var,
                                family_type = c("poisson", "nbinom2"),
                                zero_inflated = FALSE,
                                zi_formula = NULL,
                                offset_values = NULL,
                                n_quad = 40, seed = 2026) {
  family_type <- match.arg(family_type)

  set.seed(seed)
  data <- data[order(data[[id_var]], data[[order_var]]), ]

  beta_hat <- glmmTMB::fixef(fit)$cond
  X <- stats::model.matrix(fixed_formula, data)

  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]

  missing_beta <- setdiff(names(glmmTMB::fixef(fit)$cond), colnames(X))
  if (length(missing_beta) > 0) {
    stop("Missing fixed-effect columns in conditional X: ", paste(missing_beta, collapse = ", "))
  }

  if (zero_inflated) {
    beta_zi <- glmmTMB::fixef(fit)$zi
    if (is.null(zi_formula)) stop("zi_formula is required for zero-inflated RR.")
    X_zi <- stats::model.matrix(zi_formula, data)

    common_zi <- intersect(colnames(X_zi), names(beta_zi))
    X_zi <- X_zi[, common_zi, drop = FALSE]
    beta_zi <- beta_zi[common_zi]

    missing_zi <- setdiff(names(glmmTMB::fixef(fit)$zi), colnames(X_zi))
    if (length(missing_zi) > 0) {
      stop("Missing fixed-effect columns in ZI X: ", paste(missing_zi, collapse = ", "))
    }

    zi_prob_all <- as.numeric(stats::plogis(X_zi %*% beta_zi))
  } else {
    zi_prob_all <- rep(0, nrow(data))
  }

  theta <- if (family_type == "nbinom2") get_nb2_theta(fit) else NULL
  sd_b <- get_cond_re_sd(fit, id_var = id_var)

  gh <- statmod::gauss.quad.prob(n = n_quad, dist = "normal")
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights

  y <- data[[y_var]]
  ids <- data[[id_var]]

  if (is.null(offset_values)) {
    offset_values <- rep(0, nrow(data))
  } else {
    offset_values <- offset_values[order(data[[id_var]], data[[order_var]])]
  }

  rr <- numeric(nrow(data))
  u <- numeric(nrow(data))
  pred_mean <- numeric(nrow(data))
  zi_prob <- zi_prob_all

  for (id in unique(ids)) {
    rows <- which(ids == id)
    log_hist <- rep(0, length(b_nodes))

    for (pos in seq_along(rows)) {
      row <- rows[pos]

      logw <- log(base_w) + log_hist
      post_w <- normalize_log_weights(logw)

      eta_nodes <- as.numeric(X[row, , drop = FALSE] %*% beta_hat) + b_nodes + offset_values[row]
      lambda_nodes <- exp(eta_nodes)
      zi_nodes <- rep(zi_prob_all[row], length(lambda_nodes))

      k <- y[row]

      G_low <- sum(post_w * count_cdf(k - 1, lambda_nodes,
                                      family_type = family_type,
                                      theta = theta,
                                      zi = zi_nodes))

      G_high <- sum(post_w * count_cdf(k, lambda_nodes,
                                       family_type = family_type,
                                       theta = theta,
                                       zi = zi_nodes))

      a <- stats::runif(1)
      u[row] <- G_low + a * (G_high - G_low)
      rr[row] <- safe_qnorm(u[row])
      pred_mean[row] <- sum(post_w * (1 - zi_nodes) * lambda_nodes)

      pmf_nodes <- count_pmf(k, lambda_nodes,
                             family_type = family_type,
                             theta = theta,
                             zi = zi_nodes)
      pmf_nodes <- pmax(pmf_nodes, .Machine$double.xmin)
      log_hist <- log_hist + log(pmf_nodes)
    }
  }

  data.frame(data,
             rr = rr,
             u = u,
             pred_mean = pred_mean,
             zi_prob = zi_prob)
}

# ----------------------------------------------------------------------
# run_count_rr_model(): Observed recursive residuals, refit calibration and output files for
#  one count model.
# Arguments:
#   dataset_name: dataset label used in output file names
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
#   offset_expr: function of the data returning offset values, or NULL
#   index_vars: character vector of index columns for the cumulative processes
# Returns: one-row data frame summarizing the run (p_combo, p_min, leading statistic,
#  failures)
run_count_rr_model <- function(dataset_name, model_name, fit, data, fixed_formula,
                               refit_fun, id_var, order_var, y_var,
                               family_type, zero_inflated = FALSE,
                               zi_formula = NULL,
                               offset_expr = NULL,
                               index_vars) {
  cat("\n============================================================\n")
  cat("Count RR example:", dataset_name, model_name, "\n")
  cat("============================================================\n")

  offset_values <- if (is.null(offset_expr)) NULL else offset_expr(data)

  rr_obs <- rr_count_glmmTMB_ri(
    fit = fit,
    data = data,
    fixed_formula = fixed_formula,
    id_var = id_var,
    order_var = order_var,
    y_var = y_var,
    family_type = family_type,
    zero_inflated = zero_inflated,
    zi_formula = zi_formula,
    offset_values = offset_values,
    n_quad = N_QUAD,
    seed = 2026
  )

  obs_stats <- component_stats(rr_obs, index_vars, id_var, order_var)

  sim_stats <- vector("list", B_ENV)
  refit_fail <- 0

  for (b in seq_len(B_ENV)) {
    if (b %% 50 == 0) cat(dataset_name, model_name, "RR simulation", b, "of", B_ENV, "\n")

    ysim <- simulate(fit, nsim = 1, seed = 40000 + b)[[1]]
    dat_b <- data
    dat_b[[y_var]] <- as.numeric(ysim)

    fit_b <- try(suppressWarnings(refit_fun(dat_b)), silent = TRUE)

    if (inherits(fit_b, "try-error")) {
      refit_fail <- refit_fail + 1
      fit_b <- fit
    }

    offset_b <- if (is.null(offset_expr)) NULL else offset_expr(dat_b)

    rr_b <- rr_count_glmmTMB_ri(
      fit = fit_b,
      data = dat_b,
      fixed_formula = fixed_formula,
      id_var = id_var,
      order_var = order_var,
      y_var = y_var,
      family_type = family_type,
      zero_inflated = zero_inflated,
      zi_formula = zi_formula,
      offset_values = offset_b,
      n_quad = N_QUAD,
      seed = 50000 + b
    )

    sim_stats[[b]] <- component_stats(rr_b, index_vars, id_var, order_var)
  }

  sim_stats <- dplyr::bind_rows(sim_stats)
  rr_cal <- calibrate_components(obs_stats, sim_stats, ALPHA)

  prefix <- paste0(dataset_name, "_", model_name)

  readr::write_csv(rr_cal, file.path(OUT, "data", paste0(prefix,
                                                         "_RR_B1000_component_calibration.csv")))
  readr::write_csv(rr_obs, file.path(OUT, "data", paste0(prefix, "_RR_observed_residuals.csv")))
  readr::write_csv(sim_stats, file.path(OUT, "data", paste0(prefix,
                                                            "_RR_B1000_simulated_stats.csv")))

  for (v in index_vars[seq_len(min(4, length(index_vars)))]) {
    if (v %in% names(rr_obs)) {
      plot_rr_process(rr_obs, v,
                      file.path(OUT, "figures", paste0(prefix, "_RR_process_", make.names(v),
                                                       ".png")),
                      paste(prefix, "RR process by", v))
    }
  }

  summary <- data.frame(
    dataset = dataset_name,
    model_name = model_name,
    family = family_type,
    zero_inflated = zero_inflated,
    B_ENV = B_ENV,
    N_QUAD = N_QUAD,
    refit_failures = refit_fail,
    p_combo_minP = rr_cal$p_combo_minP[1],
    p_min_raw = rr_cal$p_min_raw[1],
    best_stat = rr_cal$best_stat[1],
    detected_combo = rr_cal$detected_combo[1]
  )

  readr::write_csv(summary, file.path(OUT, "data", paste0(prefix, "_RR_run_summary.csv")))
  summary
}

# ----------------------------------------------------------------------
# get_gamma_phi(): Gamma dispersion parameter of a glmmTMB fit (sigma(fit)).
# Arguments:
#   fit: fitted GLMM object
# Returns: numeric scalar
get_gamma_phi <- function(fit) {
  phi <- tryCatch(as.numeric(sigma(fit)), error = function(e) NA_real_)
  if (!is.finite(phi) || phi <= 0) {
    stop("Could not extract valid Gamma dispersion via sigma(fit).")
  }
  phi
}

# ----------------------------------------------------------------------
# rr_gamma_glmmTMB_ri(): Sequential posterior-predictive recursive residuals for a Gamma
#  random-intercept glmmTMB fit.
# Arguments:
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   n_quad: number of Gauss-Hermite quadrature nodes
#   seed: random-number seed
# Returns: data frame (sorted by cluster and order) with added columns rr, u and pred_mean
rr_gamma_glmmTMB_ri <- function(fit, data, fixed_formula,
                                id_var, order_var, y_var,
                                n_quad = 40, seed = 2026) {
  set.seed(seed)
  data <- data[order(data[[id_var]], data[[order_var]]), ]

  beta_hat <- glmmTMB::fixef(fit)$cond
  X <- stats::model.matrix(fixed_formula, data)

  common_cols <- intersect(colnames(X), names(beta_hat))
  X <- X[, common_cols, drop = FALSE]
  beta_hat <- beta_hat[common_cols]

  phi <- get_gamma_phi(fit)
  shape <- 1 / phi

  sd_b <- get_cond_re_sd(fit, id_var = id_var)

  gh <- statmod::gauss.quad.prob(n = n_quad, dist = "normal")
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights

  y <- data[[y_var]]
  ids <- data[[id_var]]

  rr <- numeric(nrow(data))
  u <- numeric(nrow(data))
  pred_mean <- numeric(nrow(data))

  for (id in unique(ids)) {
    rows <- which(ids == id)
    log_hist <- rep(0, length(b_nodes))

    for (pos in seq_along(rows)) {
      row <- rows[pos]

      logw <- log(base_w) + log_hist
      post_w <- normalize_log_weights(logw)

      eta_nodes <- as.numeric(X[row, , drop = FALSE] %*% beta_hat) + b_nodes
      mu_nodes <- exp(eta_nodes)
      scale_nodes <- mu_nodes / shape

      yy <- y[row]

      G_high <- sum(post_w * stats::pgamma(yy, shape = shape, scale = scale_nodes))
      u[row] <- pmin(pmax(G_high, .Machine$double.eps), 1 - .Machine$double.eps)
      rr[row] <- safe_qnorm(u[row])
      pred_mean[row] <- sum(post_w * mu_nodes)

      dens_nodes <- stats::dgamma(yy, shape = shape, scale = scale_nodes)
      dens_nodes <- pmax(dens_nodes, .Machine$double.xmin)
      log_hist <- log_hist + log(dens_nodes)
    }
  }

  data.frame(data, rr = rr, u = u, pred_mean = pred_mean)
}

# ----------------------------------------------------------------------
# run_gamma_rr_model(): Observed recursive residuals, refit calibration and output files for
#  one Gamma model.
# Arguments:
#   dataset_name: dataset label used in output file names
#   model_name: model label used in output file names
#   fit: fitted GLMM object
#   data: data frame used to fit the model
#   fixed_formula: one-sided formula defining the fixed-effect design matrix
#   refit_fun: function that refits the model to a data frame with a simulated response
#   id_var: name of the cluster identifier column
#   order_var: name of the within-cluster order column
#   y_var: name of the response column
#   index_vars: character vector of index columns for the cumulative processes
# Returns: one-row data frame summarizing the run
run_gamma_rr_model <- function(dataset_name, model_name, fit, data, fixed_formula,
                               refit_fun, id_var, order_var, y_var,
                               index_vars) {
  cat("\n============================================================\n")
  cat("Gamma RR example:", dataset_name, model_name, "\n")
  cat("============================================================\n")

  rr_obs <- rr_gamma_glmmTMB_ri(
    fit = fit,
    data = data,
    fixed_formula = fixed_formula,
    id_var = id_var,
    order_var = order_var,
    y_var = y_var,
    n_quad = N_QUAD,
    seed = 2026
  )

  obs_stats <- component_stats(rr_obs, index_vars, id_var, order_var)

  sim_stats <- vector("list", B_ENV)
  refit_fail <- 0

  for (b in seq_len(B_ENV)) {
    if (b %% 50 == 0) cat(dataset_name, model_name, "RR simulation", b, "of", B_ENV, "\n")

    ysim <- simulate(fit, nsim = 1, seed = 60000 + b)[[1]]
    dat_b <- data
    dat_b[[y_var]] <- as.numeric(ysim)

    fit_b <- try(suppressWarnings(refit_fun(dat_b)), silent = TRUE)

    if (inherits(fit_b, "try-error")) {
      refit_fail <- refit_fail + 1
      fit_b <- fit
    }

    rr_b <- rr_gamma_glmmTMB_ri(
      fit = fit_b,
      data = dat_b,
      fixed_formula = fixed_formula,
      id_var = id_var,
      order_var = order_var,
      y_var = y_var,
      n_quad = N_QUAD,
      seed = 70000 + b
    )

    sim_stats[[b]] <- component_stats(rr_b, index_vars, id_var, order_var)
  }

  sim_stats <- dplyr::bind_rows(sim_stats)
  rr_cal <- calibrate_components(obs_stats, sim_stats, ALPHA)

  prefix <- paste0(dataset_name, "_", model_name)

  readr::write_csv(rr_cal, file.path(OUT, "data", paste0(prefix,
                                                         "_RR_B1000_component_calibration.csv")))
  readr::write_csv(rr_obs, file.path(OUT, "data", paste0(prefix, "_RR_observed_residuals.csv")))
  readr::write_csv(sim_stats, file.path(OUT, "data", paste0(prefix,
                                                            "_RR_B1000_simulated_stats.csv")))

  for (v in index_vars[seq_len(min(4, length(index_vars)))]) {
    if (v %in% names(rr_obs)) {
      plot_rr_process(rr_obs, v,
                      file.path(OUT, "figures", paste0(prefix, "_RR_process_", make.names(v),
                                                       ".png")),
                      paste(prefix, "RR process by", v))
    }
  }

  summary <- data.frame(
    dataset = dataset_name,
    model_name = model_name,
    family = "Gamma",
    B_ENV = B_ENV,
    N_QUAD = N_QUAD,
    refit_failures = refit_fail,
    gamma_phi_assumption = "phi = sigma(fit); shape = 1/phi",
    p_combo_minP = rr_cal$p_combo_minP[1],
    p_min_raw = rr_cal$p_min_raw[1],
    best_stat = rr_cal$best_stat[1],
    detected_combo = rr_cal$detected_combo[1]
  )

  readr::write_csv(summary, file.path(OUT, "data", paste0(prefix, "_RR_run_summary.csv")))
  summary
}

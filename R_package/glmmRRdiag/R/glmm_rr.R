
#' Posterior-predictive recursive residual diagnostics for a GLMM
#'
#' Computes posterior-predictive recursive residuals for supported random-intercept
#' generalized linear mixed models and optionally performs parametric refit
#' calibration. The current manuscript prototype supports `lme4::glmer`-style
#' random-intercept models with Bernoulli/binomial, Poisson, and experimental
#' Gamma response families.
#'
#' @param fit A fitted `lme4` GLMM object inheriting from `merMod`.
#' @param data Optional data frame used to fit the model. If omitted, the model
#'   frame is extracted from `fit` when possible.
#' @param group Character name of the clustering variable. If omitted, the first
#'   random-effect grouping factor in `fit` is used.
#' @param order Optional character name of the within-cluster ordering variable.
#'   If omitted, the row order within each cluster is used.
#' @param index_vars Optional character vector of variables used for indexed
#'   cumulative residual processes. The fitted mean and within-cluster order are
#'   always included.
#' @param B Number of parametric refit-calibration replicates. Use `B = 0` to
#'   compute residuals and observed statistics only.
#' @param n_quad Number of normal quadrature nodes used for random-intercept
#'   posterior integration.
#' @param seed Optional random seed.
#' @param alpha Nominal decision level used for the `detected` field.
#' @param progress Logical; if `TRUE`, prints calibration progress.
#' @return An object of class `glmm_rr`.
#' @export
#'
#' @examples
#' \dontrun{
#' library(lme4)
#' data(cbpp, package = "lme4")
#' fit <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#'              family = binomial, data = cbpp)
#' rr <- glmm_rr(fit, group = "herd", order = "period", B = 50)
#' summary(rr)
#' plot(rr, index = "period")
#' }
glmm_rr <- function(fit,
                    data = NULL,
                    group = NULL,
                    order = NULL,
                    index_vars = NULL,
                    B = 0,
                    n_quad = 25,
                    seed = NULL,
                    alpha = 0.05,
                    progress = interactive()) {
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(fit, "merMod")) {
    stop("glmm_rr() currently supports lme4 merMod/glmer fits in this prototype.", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for glmm_rr().", call. = FALSE)
  }
  if (is.null(data)) {
    data <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
    if (is.null(data)) stop("Could not extract model frame; please supply data.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (is.null(group)) {
    fl <- lme4::getME(fit, "flist")
    if (length(fl) < 1) stop("Could not infer random-effect group; supply group=.", call. = FALSE)
    group <- names(fl)[1]
  }
  if (!group %in% names(data)) {
    # Fallback: use grouping factor from the fitted object when model.frame lacks it.
    fl <- lme4::getME(fit, "flist")
    if (group %in% names(fl)) data[[group]] <- as.factor(fl[[group]])
  }
  if (!group %in% names(data)) stop("group variable not found in data/model frame.", call. = FALSE)
  if (is.null(order)) {
    data$.rr_order <- ave(seq_len(nrow(data)), data[[group]], FUN = seq_along)
    order <- ".rr_order"
  } else if (!order %in% names(data)) {
    stop("order variable not found in data.", call. = FALSE)
  }
  y <- .rr_response(fit, data)
  fam <- stats::family(fit)
  eta_fixed <- as.numeric(stats::predict(fit, re.form = NA, type = "link"))
  mu_fixed <- fam$linkinv(eta_fixed)
  data$.rr_fitted_mean <- mu_fixed
  data$.rr_row <- seq_len(nrow(data))

  rr_data <- .compute_rr_random_intercept(
    fit = fit,
    data = data,
    y = y,
    group = group,
    order = order,
    n_quad = n_quad,
    seed = seed
  )

  idx <- unique(c(order, ".rr_fitted_mean", index_vars))
  idx <- idx[idx %in% names(rr_data)]
  observed_stats <- .component_stats(rr_data, idx, group, order)

  sim_stats <- NULL
  component_p <- NULL
  p_min <- NA_real_
  p_combo <- NA_real_
  detected <- NA
  calibration_successful <- 0L

  if (B > 0) {
    sim_list <- vector("list", B)
    used <- 0L
    response_name <- .response_name(fit)
    for (b in seq_len(B)) {
      if (progress && (b %% max(1, floor(B / 10)) == 0)) message("Calibration replicate ", b, " / ", B)
      one <- tryCatch({
        ysim <- stats::simulate(fit, nsim = 1, seed = if (is.null(seed)) NULL else seed + b)[[1]]
        dat_b <- data
        if (is.null(response_name) || !response_name %in% names(dat_b)) {
          stop("Cannot update response for this model type.")
        }
        dat_b[[response_name]] <- ysim
        fit_b <- tryCatch(stats::update(fit, data = dat_b), error = function(e) fit)
        rr_b <- .compute_rr_random_intercept(
          fit = fit_b,
          data = dat_b,
          y = .rr_response(fit_b, dat_b),
          group = group,
          order = order,
          n_quad = n_quad,
          seed = if (is.null(seed)) NULL else seed + 10000 + b
        )
        .component_stats(rr_b, idx, group, order)
      }, error = function(e) NULL)
      if (!is.null(one)) {
        used <- used + 1L
        sim_list[[used]] <- one
      }
    }
    calibration_successful <- used
    if (used > 0) {
      sim_stats <- do.call(rbind, sim_list[seq_len(used)])
      cal <- .calibrate_stats(observed_stats, sim_stats)
      component_p <- cal$component_p
      p_min <- cal$p_min
      p_combo <- cal$p_combo
      detected <- is.finite(p_combo) && p_combo < alpha
    }
  }

  out <- list(
    call = match.call(),
    fit_class = class(fit),
    family = fam$family,
    group = group,
    order = order,
    index_vars = idx,
    residuals = rr_data,
    observed_stats = observed_stats,
    sim_stats = sim_stats,
    component_p = component_p,
    p_min = p_min,
    p_combo = p_combo,
    detected = detected,
    calibration_successful = calibration_successful,
    B = B,
    n_quad = n_quad,
    alpha = alpha
  )
  class(out) <- "glmm_rr"
  out
}

#' @export
print.glmm_rr <- function(x, ...) {
  cat("Posterior-predictive recursive residual diagnostic\n")
  cat("Family:", x$family, "\n")
  cat("Group:", x$group, "\n")
  cat("Residuals:", nrow(x$residuals), "\n")
  if (x$B > 0) {
    cat("Calibration replicates requested:", x$B, "\n")
    cat("Calibration replicates used:", x$calibration_successful, "\n")
    cat("p_min:", format(x$p_min, digits = 4), "\n")
    cat("p_combo:", format(x$p_combo, digits = 4), "\n")
  } else {
    cat("No refit calibration requested (B = 0).\n")
  }
  invisible(x)
}

#' Summarize a GLMM recursive-residual diagnostic
#'
#' @param object An object returned by [glmm_rr()].
#' @param ... Ignored.
#' @return A list with observed statistics, component p-values, `p_min`, and `p_combo`.
#' @export
summary.glmm_rr <- function(object, ...) {
  tab <- data.frame(
    statistic = names(object$observed_stats),
    observed = as.numeric(object$observed_stats),
    row.names = NULL
  )
  if (!is.null(object$component_p)) {
    tab$component_p <- as.numeric(object$component_p[tab$statistic])
  }
  ans <- list(
    call = object$call,
    family = object$family,
    group = object$group,
    order = object$order,
    B = object$B,
    calibration_successful = object$calibration_successful,
    component_table = tab,
    p_min = object$p_min,
    p_combo = object$p_combo,
    detected = object$detected
  )
  class(ans) <- "summary.glmm_rr"
  ans
}

#' @export
print.summary.glmm_rr <- function(x, ...) {
  cat("Summary of posterior-predictive recursive residual diagnostic\n")
  cat("Family:", x$family, "\n")
  cat("Group:", x$group, "\n")
  cat("Calibration replicates used:", x$calibration_successful, "\n\n")
  print(x$component_table, row.names = FALSE)
  cat("\np_min:", format(x$p_min, digits = 4), "\n")
  cat("p_combo:", format(x$p_combo, digits = 4), "\n")
  cat("Detected at nominal level:", x$detected, "\n")
  invisible(x)
}

#' Plot an indexed recursive-residual process
#'
#' @param x An object returned by [glmm_rr()].
#' @param index Character name of an index in `x$index_vars`. Defaults to the
#'   first available index.
#' @param ... Additional arguments passed to `plot()`.
#' @return Invisibly returns the plotted process data.
#' @export
plot.glmm_rr <- function(x, index = NULL, ...) {
  if (is.null(index)) index <- x$index_vars[1]
  if (!index %in% names(x$residuals)) stop("Requested index not found in residual data.", call. = FALSE)
  proc <- .rr_process(x$residuals, index)
  graphics::plot(
    proc$q, proc$W, type = "l",
    xlab = index,
    ylab = expression(hat(W)[N](q)),
    main = paste("Recursive residual process indexed by", index),
    ...
  )
  graphics::abline(h = 0, lty = 2)
  invisible(proc)
}

.rr_response <- function(fit, data) {
  mf <- tryCatch(stats::model.frame(fit), error = function(e) data)
  y <- stats::model.response(mf)
  if (is.matrix(y)) stop("Matrix/binomial-count responses are not yet supported in this prototype; use Bernoulli 0/1 responses.", call. = FALSE)
  if (is.factor(y)) y <- as.numeric(y) - 1
  as.numeric(y)
}

.response_name <- function(fit) {
  f <- stats::formula(fit)
  lhs <- f[[2]]
  if (is.symbol(lhs)) return(as.character(lhs))
  NULL
}

.normal_quadrature <- function(n) {
  u <- (seq_len(n) - 0.5) / n
  list(nodes = stats::qnorm(u), weights = rep(1 / n, n))
}

.normalize_log_weights <- function(logw) {
  m <- max(logw)
  w <- exp(logw - m)
  w / sum(w)
}

.safe_qnorm <- function(u) {
  eps <- .Machine$double.eps
  stats::qnorm(pmin(pmax(u, eps), 1 - eps))
}

.compute_rr_random_intercept <- function(fit, data, y, group, order, n_quad, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fam <- stats::family(fit)
  data <- data[order(data[[group]], data[[order]], data$.rr_row), , drop = FALSE]
  y <- y[match(data$.rr_row, seq_along(y))]
  eta_fixed <- as.numeric(stats::predict(fit, newdata = data, re.form = NA, type = "link", allow.new.levels = TRUE))
  vc <- lme4::VarCorr(fit)
  sd_b <- tryCatch(as.numeric(attr(vc[[group]], "stddev"))[1], error = function(e) NA_real_)
  if (!is.finite(sd_b) || sd_b < 1e-8) sd_b <- 1e-8
  gh <- .normal_quadrature(n_quad)
  b_nodes <- sd_b * gh$nodes
  base_w <- gh$weights
  rr <- numeric(nrow(data))
  u <- numeric(nrow(data))
  pred_mean <- numeric(nrow(data))
  ids <- data[[group]]
  for (id in unique(ids)) {
    rows <- which(ids == id)
    log_hist <- rep(0, length(b_nodes))
    for (pos in seq_along(rows)) {
      row <- rows[pos]
      logw <- log(base_w) + log_hist
      post_w <- .normalize_log_weights(logw)
      eta_nodes <- eta_fixed[row] + b_nodes
      cdf <- .conditional_cdf_and_loglik(fam, y[row], eta_nodes, discrete_randomized = TRUE)
      pred_mean[row] <- sum(post_w * cdf$mean)
      G_low <- sum(post_w * cdf$G_low)
      G_high <- sum(post_w * cdf$G_high)
      a <- stats::runif(1)
      u[row] <- G_low + a * (G_high - G_low)
      rr[row] <- .safe_qnorm(u[row])
      log_hist <- log_hist + cdf$loglik
    }
  }
  data$rr <- rr
  data$u <- u
  data$rr_fitted_mean_pp <- pred_mean
  rownames(data) <- NULL
  data
}

.conditional_cdf_and_loglik <- function(fam, y, eta_nodes, discrete_randomized = TRUE) {
  family <- fam$family
  mu <- fam$linkinv(eta_nodes)
  if (family %in% c("binomial", "quasibinomial")) {
    if (!all(y %in% c(0, 1))) stop("Only Bernoulli 0/1 binomial responses are supported in this prototype.", call. = FALSE)
    p <- pmin(pmax(mu, .Machine$double.eps), 1 - .Machine$double.eps)
    if (y == 0) {
      G_low <- rep(0, length(p)); G_high <- 1 - p; loglik <- stats::dbinom(0, 1, p, log = TRUE)
    } else {
      G_low <- 1 - p; G_high <- rep(1, length(p)); loglik <- stats::dbinom(1, 1, p, log = TRUE)
    }
    return(list(G_low = G_low, G_high = G_high, loglik = loglik, mean = p))
  }
  if (family %in% c("poisson", "quasipoisson")) {
    lambda <- pmax(mu, .Machine$double.eps)
    yy <- as.integer(round(y))
    G_low <- if (yy <= 0) rep(0, length(lambda)) else stats::ppois(yy - 1, lambda)
    G_high <- stats::ppois(yy, lambda)
    loglik <- stats::dpois(yy, lambda, log = TRUE)
    return(list(G_low = G_low, G_high = G_high, loglik = loglik, mean = lambda))
  }
  if (family == "Gamma") {
    # Experimental continuous Gamma residual. Dispersion extraction is deliberately
    # conservative because GLMM software differs in how it reports scale.
    dispersion <- 1
    shape <- 1 / dispersion
    scale <- pmax(mu / shape, .Machine$double.eps)
    yy <- pmax(y, .Machine$double.eps)
    G <- stats::pgamma(yy, shape = shape, scale = scale)
    loglik <- stats::dgamma(yy, shape = shape, scale = scale, log = TRUE)
    return(list(G_low = G, G_high = G, loglik = loglik, mean = mu))
  }
  stop("Family not supported in this prototype: ", family, call. = FALSE)
}

.rr_process <- function(rr_data, q_var) {
  ord <- order(rr_data[[q_var]], rr_data$.rr_row)
  q <- rr_data[[q_var]][ord]
  r <- rr_data$rr[ord]
  data.frame(q = q, W = cumsum(r) / sqrt(length(r)))
}

.component_stats <- function(rr_data, index_vars, group, order) {
  out <- list()
  for (v in index_vars) {
    if (!v %in% names(rr_data)) next
    if (length(unique(rr_data[[v]])) < 2) next
    proc <- .rr_process(rr_data, v)
    nm <- make.names(v)
    out[[paste0("Tinf_", nm)]] <- max(abs(proc$W), na.rm = TRUE)
    out[[paste0("T2_", nm)]] <- mean(proc$W^2, na.rm = TRUE)
  }
  r <- rr_data$rr
  out$Tmax_r <- max(abs(r), na.rm = TRUE)
  C <- stats::aggregate(rr_data$rr, list(rr_data[[group]]), function(z) sum(z) / sqrt(length(z)))
  out$Tmax_C <- max(abs(C$x), na.rm = TRUE)
  out$Tdisp <- abs(mean(r^2, na.rm = TRUE) - 1)
  lagdat <- rr_data[order(rr_data[[group]], rr_data[[order]], rr_data$.rr_row), , drop = FALSE]
  lagdat$lag_rr <- ave(lagdat$rr, lagdat[[group]], FUN = function(z) c(NA, z[-length(z)]))
  lagdat <- lagdat[!is.na(lagdat$lag_rr), , drop = FALSE]
  out$Tlag <- if (nrow(lagdat) > 2) abs(stats::cor(lagdat$rr, lagdat$lag_rr)) else NA_real_
  unlist(out)
}

.p_from_sim <- function(obs, sim) {
  valid <- sim[is.finite(sim)]
  if (!is.finite(obs) || length(valid) == 0) return(NA_real_)
  (1 + sum(valid >= obs)) / (length(valid) + 1)
}

.calibrate_stats <- function(obs, sim_stats) {
  names_common <- intersect(names(obs), colnames(sim_stats))
  names_common <- names_common[is.finite(obs[names_common])]
  component_p <- stats::setNames(rep(NA_real_, length(names_common)), names_common)
  for (s in names_common) component_p[s] <- .p_from_sim(obs[s], sim_stats[[s]])
  p_min <- min(component_p, na.rm = TRUE)
  sim_p <- matrix(NA_real_, nrow = nrow(sim_stats), ncol = length(names_common))
  colnames(sim_p) <- names_common
  for (j in seq_along(names_common)) {
    vals <- sim_stats[[names_common[j]]]
    sim_p[, j] <- vapply(vals, .p_from_sim, numeric(1), sim = vals)
  }
  sim_min <- apply(sim_p, 1, min, na.rm = TRUE)
  p_combo <- (1 + sum(sim_min <= p_min, na.rm = TRUE)) / (length(sim_min) + 1)
  list(component_p = component_p, p_min = p_min, p_combo = p_combo)
}

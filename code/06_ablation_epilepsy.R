# ============================================================
# 06_ablation_epilepsy.R
#
# One-factor ablation of the recursive-residual (RR) diagnostic on the
# MASS::epil data (response to Reviewer 1, second round).
#
# Part A  Worked example. The three epilepsy models of Section 6.2
#         (Poisson RI, NB2 RI, refined NB2) are diagnosed with one
#         building block replaced at a time; the other two are held at
#         the proposed choice.
#           (i)   residual construction: sequential posterior (proposed)
#                 vs marginal (prior of b_i, no updating) vs posterior
#                 given the full cluster response vector;
#           (ii)  reference distribution: refit (proposed) vs plug-in
#                 (simulate at theta-hat, no refit) vs known-parameter
#                 (i.i.d. N(0,1) residuals; Propositions 3.1-3.5);
#           (iii) statistics: full component set (proposed) vs indexed
#                 cumulative statistics only vs direct distributional
#                 statistics only (KS, Cramer-von Mises, dispersion,
#                 max |r|), all refit-calibrated, and the classical
#                 uncalibrated Kolmogorov-Smirnov and Shapiro-Wilk tests.
#         All variants share the fitted models, the simulated datasets
#         and the PIT randomization draws (common random numbers).
#         The "proposed" row reuses the seeds of scripts 03 and 05 and
#         should reproduce the published p_combo values (0.001, 0.001,
#         0.166) up to package-version differences.
#
# Part B  Design-anchored replication. The same variants are applied to
#         datasets generated from the fitted epilepsy models on the
#         observed design, giving rejection rates under a null, a
#         mean-structure misspecification, and serial dependence.
#
# Run from the root of the code/data supplement:
#   source("code/06_ablation_epilepsy.R")
# Recommended: first run with SMOKE_TEST <- TRUE (about 1-3 minutes).
# Outputs: ./results/ablation_epilepsy/{data,tables,figures,logs}
# ============================================================

required <- c("glmmTMB", "statmod", "MASS")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Please install missing packages first: ",
                                   paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages(library(glmmTMB))

# ---------------- Settings ----------------
SMOKE_TEST <- FALSE   # TRUE: tiny end-to-end check before the full run
RUN_PART_A <- TRUE    # worked example on the observed data
RUN_PART_B <- TRUE    # design-anchored rejection rates

N_QUAD     <- 40
ALPHA      <- 0.05
B_A        <- if (SMOKE_TEST) 19L else 1000L   # reference replicates, Part A
N_REP      <- if (SMOKE_TEST) 3L  else 200L    # generated datasets per scenario, Part B
B_INNER    <- if (SMOKE_TEST) 19L else 199L    # reference replicates per dataset, Part B
RHO_SERIAL <- c(0.25, 0.50)                    # lagged-count coefficients, Part B
N_CORES <- as.integer(2)

OUT <- file.path(getwd(), "results",
                 if (SMOKE_TEST) "ablation_epilepsy_SMOKE" else "ablation_epilepsy")
for (sub in c("data", "tables", "figures", "logs")) dir.create(file.path(OUT, sub),
    recursive = TRUE, showWarnings = FALSE)
cat("Outputs will be written to:", OUT, "\n")
cat("Part B cores:", N_CORES, "\n")

# ---------------- Data (identical preparation to scripts 03 and 05) ----------------
data(epil, package = "MASS")
ep <- epil
ep$subject    <- factor(ep$subject)
ep$period_num <- as.numeric(as.character(ep$period))
ep$period_fac <- factor(ep$period)
ep$trt        <- factor(ep$trt)
ep$log_base   <- log(ep$base + 1)
ep$log_base_s <- as.numeric(scale(ep$log_base))
ep$age_s      <- as.numeric(scale(ep$age))

SPECS <- list(
  Poisson_RI = list(
    label = "Poisson RI", family = "poisson",
    formula = y ~ trt + base + age + period_fac + (1 | subject),
    fixed = ~ trt + base + age + period_fac,
    index = c("period_num", "pred_mean", "base", "age", "trt"),
    obs_seed = 2026, sim_seed = 40000, res_seed = 50000, kp_seed = 700000),
  NB2_RI = list(
    label = "NB2 RI", family = "nbinom2",
    formula = y ~ trt + base + age + period_fac + (1 | subject),
    fixed = ~ trt + base + age + period_fac,
    index = c("period_num", "pred_mean", "base", "age", "trt"),
    obs_seed = 2026, sim_seed = 40000, res_seed = 50000, kp_seed = 700000),
  NB2_refined = list(
    label = "NB2 refined", family = "nbinom2",
    formula = y ~ trt * period_fac + poly(log_base_s, 2) + age_s + (1 | subject),
    fixed = ~ trt * period_fac + poly(log_base_s, 2) + age_s,
    index = c("period_num", "pred_mean", "base", "log_base", "age", "trt"),
    obs_seed = 2026, sim_seed = 90000, res_seed = 100000, kp_seed = 700000)
)

CONSTR <- c("seq", "marg", "full")

VARIANTS <- data.frame(
  key = c("proposed", "marginal", "fullpost", "plugin", "known", "neither",
          "indexed", "direct", "ks", "sw"),
  block = c("Proposed", "(i) Residual", "(i) Residual", "(ii) Reference", "(ii) Reference",
            "(i)+(ii) Both", "(iii) Statistics", "(iii) Statistics", "(iii) Statistics",
                "(iii) Statistics"),
  label_plain = c("Sequential residuals, refit reference, full component set",
                  "Marginal residuals (prior of b_i, no updating)",
                  "Full-cluster posterior residuals",
                  "Plug-in reference (no refit)",
                  "Known-parameter reference (iid N(0,1))",
                  "Marginal residuals with plug-in reference (both removed)",
                  "Indexed cumulative statistics only",
                  "Direct distributional statistics only",
                  "Classical Kolmogorov-Smirnov, uncalibrated",
                  "Classical Shapiro-Wilk, uncalibrated"),
  label_tex = c("Sequential residuals, refit reference, full set",
                "Marginal residuals (prior of $b_i$, no updating)",
                "Full-cluster posterior residuals",
                "Plug-in reference (simulate at $\\hat\\theta$, no refit)",
                "Known-parameter reference (i.i.d.\\ $N(0,1)$)",
                "Marginal residuals with plug-in reference",
                "Indexed cumulative statistics only",
                "Direct distributional statistics only",
                "Classical Kolmogorov--Smirnov, uncalibrated",
                "Classical Shapiro--Wilk, uncalibrated"),
  stringsAsFactors = FALSE)

# ---------------- Model fitting and parameter extraction ----------------
# Function definitions used in this script
source(file.path("code", "functions", "functions_06_ablation_epilepsy.R"), local = TRUE)


# ---------------- Residual constructions ----------------

# ---------------- Statistics ----------------


# ---------------- Calibration (identical rules to scripts 03 and 05) ----------------


# ---------------- One complete diagnosis with all variants ----------------


# ---------------- LaTeX helpers ----------------


# ============================================================
# Part A: worked example on the observed data
# ============================================================
if (RUN_PART_A) {
  resA <- list()
  for (m in names(SPECS)) {
    spec <- SPECS[[m]]
    cat("\nPart A:", spec$label, "(B =", B_A, ")\n")
    t0 <- Sys.time()
    fit <- fit_model(spec, ep)
    resA[[m]] <- diagnose(fit, ep, spec, B_A,
                          seeds = list(obs = spec$obs_seed, sim = spec$sim_seed,
                                       res = spec$res_seed, kp = spec$kp_seed),
                          keep_process = if (m == "NB2_RI") c("base",
                              "period_num") else character(0))
    cat("  finished in", format(round(Sys.time() - t0, 1)),
        "| refit failures:", resA[[m]]$n_refit_fail, "| skipped:", resA[[m]]$n_skip, "\n")
  }
  saveRDS(resA, file.path(OUT, "data", "partA_results.rds"))

  sumA <- do.call(rbind, lapply(names(resA), function(m) cbind(model = m, summ_row(resA[[m]]),
                                                               stringsAsFactors = FALSE)))
  sumA$variant <- VARIANTS$label_plain[match(sumA$key, VARIANTS$key)]
  write.csv(sumA, file.path(OUT, "data", "partA_pvalues_by_variant.csv"), row.names = FALSE)

  compA <- do.call(rbind, lapply(names(resA), function(m) do.call(rbind,
      lapply(names(resA[[m]]$V), function(k) {
    x <- resA[[m]]$V[[k]]
    data.frame(model = m, variant = k, statistic = names(x$p), p = unname(x$p),
               stringsAsFactors = FALSE)
  }))))
  write.csv(compA, file.path(OUT, "data", "partA_component_pvalues.csv"), row.names = FALSE)

  obsA <- do.call(rbind, lapply(names(resA), function(m) do.call(rbind, lapply(CONSTR, function(w) {
    o <- resA[[m]]$obs[[w]]
    data.frame(model = m, construction = w, statistic = names(o), observed = unname(o),
               stringsAsFactors = FALSE)
  }))))
  write.csv(obsA, file.path(OUT, "data", "partA_observed_statistics.csv"), row.names = FALSE)

  # Reproduction check against the published proposed p_combo values
  pub <- c(Poisson_RI = 0.000999, NB2_RI = 0.000999, NB2_refined = 0.165834)
  chk <- data.frame(model = names(pub), published = unname(pub),
                    reproduced = vapply(names(pub), function(m) resA[[m]]$V$proposed$p_combo,
                                        numeric(1)))
  write.csv(chk, file.path(OUT, "data", "partA_reproduction_check.csv"), row.names = FALSE)
  cat("\nReproduction check (proposed variant):\n"); print(chk, row.names = FALSE)

  # Table B1: p_combo by variant and model
  mods <- names(SPECS)
  L <- c("\\begin{table}[htbp]", "\\centering", "\\small",
         paste0("\\caption{One-factor ablation of the recursive-residual diagnostic on the ",
                "epilepsy data. ",
                "Each row replaces one building block and holds the other two at the ",
                "proposed choice. ",
                "Entries are calibrated omnibus $p$-values based on $B=", B_A,
                    "$ simulated reference datasets; ",
                "the last two rows are classical uncalibrated $p$-values computed on the ",
                "sequential residuals. ",
                "The last column gives the leading component for the NB2 random-intercept model.}"),
         "\\label{tab:ablation-epil}",
         "\\resizebox{\\textwidth}{!}{%",
         "\\begin{tabular}{llcccl}", "\\toprule",
         "Block & Variant & Poisson RI & NB2 RI & NB2 refined & Leading (NB2 RI) \\\\", "\\midrule")
  blk <- block_col(VARIANTS$block)
  for (i in seq_len(nrow(VARIANTS))) {
    k <- VARIANTS$key[i]
    ps <- vapply(mods, function(m) fmt3(sumA$p[sumA$model == m & sumA$key == k]), character(1))
    ld <- sumA$lead[sumA$model == "NB2_RI" & sumA$key == k]
    ld <- if (length(ld) == 1 && !is.na(ld)) stat_tex(ld) else "--"
    L <- c(L, paste(blk[i], VARIANTS$label_tex[i], ps[1], ps[2], ps[3], paste0(ld, " \\\\"),
                    sep = " & "))
    if (k %in% c("proposed", "fullpost", "known", "neither")) L <- c(L, "\\midrule")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(OUT, "tables", "table_ablation_epil_pcombo.tex"))

  # Table B2: component p-values, NB2 RI
  V2 <- resA$NB2_RI$V
  keys2 <- c("proposed", "marginal", "fullpost", "plugin", "known")
  stats2 <- names(V2$proposed$p)
  L <- c("\\begin{table}[htbp]", "\\centering", "\\small",
         paste0("\\caption{Component-level calibrated $p$-values for the NB2 ",
                "random-intercept epilepsy model ",
                "under the three residual constructions (refit reference) and the three ",
                "reference distributions ",
                "(sequential residuals). $B=", B_A, "$.}"),
         "\\label{tab:ablation-epil-components}",
         "\\begin{tabular}{lccccc}", "\\toprule",
         paste0(" & \\multicolumn{3}{c}{(i) Residual construction} & \\multicolumn{2}{c}{(ii) ",
                "Reference} \\\\"),
         "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
         "Statistic & Sequential & Marginal & Full posterior & Plug-in & Known-par. \\\\",
             "\\midrule")
  for (s in stats2) {
    vals <- vapply(keys2,
                   function(k) fmt3(if (s %in% names(V2[[k]]$p)) V2[[k]]$p[[s]] else NA_real_),
                       character(1))
    L <- c(L, paste0(stat_tex(s), " & ", paste(vals, collapse = " & "), " \\\\"))
  }
  L <- c(L, "\\midrule",
         paste0("$p_{\\mathrm{combo}}$ & ", paste(vapply(keys2,
                                                         function(k) fmt3(V2[[k]]$p_combo),
                                                             character(1)), collapse = " & "),
                                                                 " \\\\"),
         "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(L, file.path(OUT, "tables", "table_ablation_epil_components.tex"))

  # Figure B1: indexed processes, NB2 RI, three residual constructions
  res <- resA$NB2_RI
  titles <- c(seq = "Sequential (proposed)", marg = "Marginal", full = "Full-cluster posterior")
  vkey <- c(seq = "proposed", marg = "marginal", full = "fullpost")
  grDevices::pdf(file.path(OUT, "figures", "fig_ablation_epil_processes.pdf"), width = 10,
                 height = 6.4)
  op <- graphics::par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.4, 0.8))
  for (v in c("base", "period_num")) for (w in CONSTR) {
    Pm <- res$P[[paste(w, v, sep = ":")]]
    Pm <- Pm[stats::complete.cases(Pm), , drop = FALSE]
    Wobs <- proc_W(res$d_obs, w, v)
    x <- seq_along(Wobs) / length(Wobs)
    lo <- apply(Pm, 2, stats::quantile, probs = 0.025)
    hi <- apply(Pm, 2, stats::quantile, probs = 0.975)
    graphics::plot(x, Wobs, type = "n", ylim = range(c(lo, hi, Wobs)),
                   xlab = paste("Proportion of observations ordered by",
                                if (v == "base") "baseline count" else "period"),
                   ylab = "Cumulative residual process", main = titles[[w]])
    graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), col = "grey85", border = NA)
    graphics::abline(h = 0, lty = 2)
    graphics::lines(x, Wobs, lwd = 1.5)
    pc <- res$V[[vkey[[w]]]]$p[[paste0("T2_", v)]]
    graphics::legend("topleft", bty = "n", cex = 0.9,
                     legend = paste0("T2 component p = ", sprintf("%.3f", pc)))
  }
  graphics::par(op)
  grDevices::dev.off()
}

# ============================================================
# Part B: design-anchored replication (rejection rates)
# ============================================================
if (RUN_PART_B) {
  fit_ref_obs <- fit_model(SPECS$NB2_refined, ep)
  fit_nb2_obs <- fit_model(SPECS$NB2_RI, ep)
  c0 <- mean(log1p(ep$y))

  # NB2 RI with a lagged-count term: log mu_it = eta_it + b_i + rho * {log(1 + y_i,t-1) -
  # c0}, t >= 2
  gen_serial <- function(rho, seed) {
    beta <- fixef(fit_nb2_obs)$cond
    X <- stats::model.matrix(SPECS$NB2_RI$fixed, ep)
    eta <- as.numeric(X[, names(beta), drop = FALSE] %*% beta)
    sdb <- re_sd(fit_nb2_obs); th <- nb2_theta(fit_nb2_obs)
    set.seed(seed)
    b <- stats::rnorm(nlevels(ep$subject), 0, sdb); names(b) <- levels(ep$subject)
    y <- numeric(nrow(ep))
    prev_y <- NA_real_; prev_id <- ""
    for (row in order(ep$subject, ep$period_num)) {
      id <- as.character(ep$subject[row])
      lp <- eta[row] + b[[id]]
      if (identical(prev_id, id)) lp <- lp + rho * (log1p(prev_y) - c0)
      y[row] <- stats::rnbinom(1, mu = exp(lp), size = th)
      prev_y <- y[row]; prev_id <- id
    }
    y
  }

  SCEN <- list(
    null = list(label = "Null: refined NB2 generated and fitted", analysis = "NB2_refined",
                gen = function(seed) as.numeric(stats::simulate(fit_ref_obs, nsim = 1,
                    seed = seed)[[1]])),
    mean = list(label = "Mean misspecification: refined NB2 generated, NB2 RI fitted",
                analysis = "NB2_RI",
                gen = function(seed) as.numeric(stats::simulate(fit_ref_obs, nsim = 1,
                    seed = seed)[[1]])),
    serial1 = list(label = paste0("Serial dependence, rho = ", RHO_SERIAL[1],
                                  ": NB2 RI fitted"), analysis = "NB2_RI",
                   gen = function(seed) gen_serial(RHO_SERIAL[1], seed)),
    serial2 = list(label = paste0("Serial dependence, rho = ", RHO_SERIAL[2],
                                  ": NB2 RI fitted"), analysis = "NB2_RI",
                   gen = function(seed) gen_serial(RHO_SERIAL[2], seed)))

  one_rep <- function(r, s_id, sc) {
    base_seed <- 1e7 * s_id + 1e4 * r
    d <- ep
    d$y <- sc$gen(base_seed + 9500)
    spec <- SPECS[[sc$analysis]]
    fit <- tryCatch(fit_model(spec, d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    res <- tryCatch(diagnose(fit, d, spec, B_INNER,
                             seeds = list(obs = base_seed + 9000, sim = base_seed,
                                          res = base_seed + 5000, kp = base_seed + 7000)),
                    error = function(e) NULL)
    if (is.null(res)) return(NULL)
    out <- summ_row(res)
    out$rep <- r
    out
  }

  resB <- list()
  for (s_id in seq_along(SCEN)) {
    nm <- names(SCEN)[s_id]; sc <- SCEN[[s_id]]
    ck <- file.path(OUT, "data", paste0("partB_", nm, ".rds"))
    if (file.exists(ck)) { cat("\nPart B:", sc$label,
                               "- loaded checkpoint\n"); resB[[nm]] <- readRDS(ck); next }
    cat("\nPart B:", sc$label, "(", N_REP, "datasets, B =", B_INNER, ")\n")
    t0 <- Sys.time()
    reps <- if (N_CORES > 1) {
      parallel::mclapply(seq_len(N_REP), one_rep, s_id = s_id, sc = sc,
                         mc.cores = N_CORES, mc.preschedule = FALSE)
    } else {
      lapply(seq_len(N_REP), one_rep, s_id = s_id, sc = sc)
    }
    reps <- reps[vapply(reps, is.data.frame, logical(1))]
    if (length(reps) == 0) stop("Part B scenario ", nm, ": no successful replicates.")
    R <- do.call(rbind, reps)
    R$scenario <- nm
    saveRDS(R, ck)
    resB[[nm]] <- R
    cat("  finished in", format(round(Sys.time() - t0, 1)), "| successful datasets:",
        length(reps), "\n")
  }

  wilson <- function(x, n, z = stats::qnorm(0.975)) {
    if (n == 0) return(c(NA_real_, NA_real_))
    p <- x / n; den <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / den
    hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
    c(ctr - hw, ctr + hw)
  }

  sumB <- do.call(rbind, lapply(names(resB), function(nm) {
    R <- resB[[nm]]
    do.call(rbind, lapply(VARIANTS$key, function(k) {
      pk <- R$p[R$key == k]; pk <- pk[is.finite(pk)]
      x <- sum(pk <= ALPHA); n <- length(pk); ci <- wilson(x, n)
      data.frame(scenario = nm, key = k, variant = VARIANTS$label_plain[VARIANTS$key == k],
                 n = n, rejections = x, rate = if (n > 0) x / n else NA_real_,
                 ci_low = ci[1], ci_high = ci[2], stringsAsFactors = FALSE)
    }))
  }))
  write.csv(sumB, file.path(OUT, "data", "partB_rejection_rates.csv"), row.names = FALSE)

  leadB <- do.call(rbind, lapply(names(resB), function(nm) {
    R <- resB[[nm]]
    Ld <- R$lead[R$key == "proposed" & is.finite(R$p) & R$p <= ALPHA]
    if (length(Ld) == 0) return(NULL)
    tb <- sort(table(Ld), decreasing = TRUE)
    data.frame(scenario = nm, leading_component = names(tb), count = as.integer(tb),
               n_rejected = length(Ld), stringsAsFactors = FALSE)
  }))
  if (!is.null(leadB)) write.csv(leadB, file.path(OUT, "data",
                                                  "partB_leading_components_proposed.csv"),
                                                      row.names = FALSE)

  # Table B3: rejection rates
  sc_names <- names(SCEN)
  heads <- c("Null", "Mean misspec.", paste0("Serial, $\\rho=", RHO_SERIAL, "$"))
  L <- c("\\begin{table}[htbp]", "\\centering", "\\small",
         paste0("\\caption{Rejection rates (\\%) at the 5\\% level for datasets generated ",
                "from the fitted epilepsy ",
                "models on the observed design, with Wilson 95\\% intervals. Null: refined ",
                "NB2 generated and fitted. ",
                "Mean misspecification: refined NB2 generated, NB2 random-intercept model fitted. ",
                "Serial dependence: NB2 random-intercept model with a lagged log-count term ",
                "generated, ",
                "NB2 random-intercept model fitted. ", N_REP, " datasets per scenario; $B=",
                    B_INNER,
                "$ reference replicates per dataset.}"),
         "\\label{tab:ablation-epil-rates}",
         "\\resizebox{\\textwidth}{!}{%",
         "\\begin{tabular}{llcccc}", "\\toprule",
         paste0("Block & Variant & ", paste(heads, collapse = " & "), " \\\\"), "\\midrule")
  blk <- block_col(VARIANTS$block)
  for (i in seq_len(nrow(VARIANTS))) {
    k <- VARIANTS$key[i]
    cells <- vapply(sc_names, function(nm) {
      z <- sumB[sumB$scenario == nm & sumB$key == k, ]
      if (nrow(z) == 0 || !is.finite(z$rate)) return("--")
      sprintf("%.1f (%.1f, %.1f)", 100 * z$rate, 100 * z$ci_low, 100 * z$ci_high)
    }, character(1))
    L <- c(L, paste0(blk[i], " & ", VARIANTS$label_tex[i], " & ", paste(cells,
        collapse = " & "), " \\\\"))
    if (k %in% c("proposed", "fullpost", "known", "neither")) L <- c(L, "\\midrule")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(L, file.path(OUT, "tables", "table_ablation_epil_rates.tex"))
}

sink(file.path(OUT, "logs", "sessionInfo_ablation_epilepsy.txt"))
print(sessionInfo())
sink()
cat("\nDONE. Outputs are in:\n", OUT, "\n")

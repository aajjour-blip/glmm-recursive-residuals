# ============================================================
# RR_GLMM manuscript -- Bernoulli Q1-level power and null-level simulation
# Self-contained R version
# ============================================================
#
# What this script does
#   1. Loads required CRAN packages.
#   2. Runs the Bernoulli random-intercept GLMM RR simulation.
#   3. Uses refit calibration for the minimum-p omnibus statistic.
#   4. Produces CSV results, LaTeX tables, figures, and logs.
#   5. Supports resume from partial results.
#
# Recommended R use
#   - Run from the root of the code/data supplement:
#       source("code/01_simulation_bernoulli_power.R")
#
# FULL Q1 settings are the default. This is computationally intensive.
# For a quick test before launching the full run, run this first:
#       Sys.setenv(RR_GLMM_RUN_MODE = "SMOKE_TEST")
#       source("code/01_simulation_bernoulli_power.R")
# Then, for the full run:
#       Sys.setenv(RR_GLMM_RUN_MODE = "FULL_Q1")
#       source("code/01_simulation_bernoulli_power.R")
#
# Output location
#   - By default, outputs are written under ./results/.
#   - To force a specific relative output root, run before sourcing:
#       Sys.setenv(RR_GLMM_OUTPUT_ROOT = "results")
#
# ============================================================

# -----------------------------
# 0. Basic session safeguards
# -----------------------------

options(stringsAsFactors = FALSE)
options(warn = 1)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Make dplyr/tidyr-style code more stable across sessions.
Sys.setenv(R_TESTS = "")

run_start_time <- Sys.time()

# -----------------------------
# 1. Load packages
# -----------------------------

required_packages <- c(
  "lme4",
  "statmod",
  "dplyr",
  "ggplot2",
  "readr",
  "tibble"
)

# Function definitions used in this script
source(file.path("code", "functions", "functions_01_simulation_bernoulli_power.R"), local = TRUE)

# Required packages are loaded below; install missing packages manually if needed.
suppressPackageStartupMessages({
  library(lme4)
  library(statmod)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

set.seed(2026)

# -----------------------------
# 2. Run mode and simulation settings
# -----------------------------

RUN_MODE <- toupper(Sys.getenv("RR_GLMM_RUN_MODE", unset = "FULL_Q1"))
if (!RUN_MODE %in% c("FULL_Q1", "SMOKE_TEST")) {
  stop("RR_GLMM_RUN_MODE must be either 'FULL_Q1' or 'SMOKE_TEST'.")
}

if (RUN_MODE == "FULL_Q1") {
  N_ALT  <- 200     # non-null datasets per design cell
  N_NULL <- 500     # null datasets per design cell, gamma = 0
  B_CAL  <- 200     # refit-calibration simulations per fitted dataset
  M_VALUES   <- c(50, 100, 200)
  N_VALUES   <- c(4, 8, 12)
  GAMMA_ALT  <- c(0.25, 0.50, 0.75, 1.00)
} else {
  # Small settings only to confirm that the script, output writing,
  # and plotting work correctly.
  N_ALT  <- 2
  N_NULL <- 2
  B_CAL  <- 50  # >= 40 required by the calibration validity gate
  M_VALUES   <- c(20)
  N_VALUES   <- c(4)
  GAMMA_ALT  <- c(0.50)
}

ALPHA <- 0.05
N_QUAD <- 25
GAMMA_NULL <- 0.00

BETA0 <- -1.00
BETA1 <-  0.80
SIGMA_B <- 0.80

RESUME <- TRUE

# -----------------------------
# 3. Output folder setup
# -----------------------------


OUTPUT_ROOT <- get_output_root()
OUT <- file.path(OUTPUT_ROOT, "bernoulli_power_results")

subdirs <- c("data", "figures", "tables_latex", "logs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (s in subdirs) dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT, "logs", "run_progress.log")


log_msg("Run mode: ", RUN_MODE)
log_msg("Output folder: ", OUT)

# -----------------------------
# 4. Utility functions
# -----------------------------


# -----------------------------
# 5. Data generation and RR calculation
# -----------------------------


# -----------------------------
# 6. Build design grids
# -----------------------------

null_grid <- expand.grid(
  m = M_VALUES,
  n_i = N_VALUES,
  gamma = GAMMA_NULL,
  rep_id = seq_len(N_NULL),
  KEEP.OUT.ATTRS = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(row_type = "null_level")

alt_grid <- expand.grid(
  m = M_VALUES,
  n_i = N_VALUES,
  gamma = GAMMA_ALT,
  rep_id = seq_len(N_ALT),
  KEEP.OUT.ATTRS = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(row_type = "power")

design_grid <- dplyr::bind_rows(null_grid, alt_grid) |>
  dplyr::arrange(row_type, gamma, m, n_i, rep_id)

design_grid$row_id <- seq_len(nrow(design_grid))

readr::write_csv(design_grid, file.path(OUT, "data", "bernoulli_Q1_design_grid.csv"))

config <- tibble::tibble(
  run_type = "RR_GLMM_manuscript_Bernoulli_Q1_power_level",
  run_mode = RUN_MODE,
  N_ALT = N_ALT,
  N_NULL = N_NULL,
  B_CAL = B_CAL,
  ALPHA = ALPHA,
  N_QUAD = N_QUAD,
  M_VALUES = paste(M_VALUES, collapse = ","),
  N_VALUES = paste(N_VALUES, collapse = ","),
  GAMMA_ALT = paste(GAMMA_ALT, collapse = ","),
  GAMMA_NULL = GAMMA_NULL,
  output_folder = OUT,
  run_start_time = as.character(run_start_time)
)

readr::write_csv(config, file.path(OUT, "logs", "run_config.csv"))

# -----------------------------
# 7. Resume and run
# -----------------------------

partial_file <- file.path(OUT, "data", "bernoulli_Q1_replicate_results_partial.csv")

completed <- tibble::tibble()
if (RESUME && file.exists(partial_file)) {
  completed <- tryCatch(readr::read_csv(partial_file, show_col_types = FALSE),
                        error = function(e) tibble::tibble())
  if (nrow(completed) > 0) {
    completed <- completed |>
      dplyr::distinct(row_type, m, n_i, gamma, rep_id, .keep_all = TRUE)
  }
}

if (nrow(completed) > 0) {
  remaining_grid <- design_grid |>
    dplyr::anti_join(
      completed |> dplyr::select(row_type, m, n_i, gamma, rep_id),
      by = c("row_type", "m", "n_i", "gamma", "rep_id")
    )
} else {
  remaining_grid <- design_grid
}

log_msg("Total design rows: ", nrow(design_grid))
log_msg("Completed rows: ", nrow(completed))
log_msg("Remaining rows: ", nrow(remaining_grid))
log_msg("Total calibration refits from scratch: ", nrow(design_grid) * B_CAL)

new_results <- vector("list", nrow(remaining_grid))

if (nrow(remaining_grid) > 0) {
  for (ii in seq_len(nrow(remaining_grid))) {
    row <- remaining_grid[ii, ]
    elapsed <- difftime(Sys.time(), run_start_time, units = "mins")

    log_msg(
      "Remaining row ", ii, " of ", nrow(remaining_grid),
      " | original row_id = ", row$row_id,
      " | type = ", row$row_type,
      " | gamma = ", row$gamma,
      " | m = ", row$m,
      " | n_i = ", row$n_i,
      " | rep = ", row$rep_id,
      " | elapsed min = ", round(as.numeric(elapsed), 2)
    )

    res <- run_one_replicate(
      m = row$m,
      n_i = row$n_i,
      gamma = row$gamma,
      rep_id = row$rep_id,
      row_type = row$row_type
    )

    new_results[[ii]] <- res

    if (!file.exists(partial_file) || file.info(partial_file)$size == 0) {
      readr::write_csv(res, partial_file)
    } else {
      readr::write_csv(res, partial_file, append = TRUE)
    }
  }
}

new_results_df <- dplyr::bind_rows(new_results)

replicate_results <- dplyr::bind_rows(completed, new_results_df) |>
  dplyr::distinct(row_type, m, n_i, gamma, rep_id, .keep_all = TRUE) |>
  dplyr::arrange(row_type, gamma, m, n_i, rep_id)

readr::write_csv(replicate_results, file.path(OUT, "data", "bernoulli_Q1_replicate_results.csv"))

# -----------------------------
# 8. Summaries with MC SE and Wilson CI
# -----------------------------

cell_summary <- replicate_results |>
  dplyr::group_by(row_type, m, n_i, gamma) |>
  dplyr::summarise(
    N_requested = ifelse(dplyr::first(row_type) == "null_level", N_NULL, N_ALT),
    N_success = sum(is.finite(p_combo_minP)),
    detected_n = sum(detected_combo, na.rm = TRUE),
    detection_rate = ifelse(N_success > 0, detected_n / N_success, NA_real_),
    mc_se = mc_se(detection_rate, N_success),
    mean_calibration_successful = mean(calibration_successful, na.rm = TRUE),
    fit_failures = sum(fit_failed, na.rm = TRUE),
    mean_refit_failures = mean(refit_failures, na.rm = TRUE),
    mean_rr_failures = mean(rr_failures, na.rm = TRUE),
    mean_singular_refits = mean(singular_refits, na.rm = TRUE),
    median_p_combo = median(p_combo_minP, na.rm = TRUE),
    most_common_best_stat = mode_string(best_stat),
    .groups = "drop"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ci_low = wilson_ci(detected_n, N_success)[1],
    ci_high = wilson_ci(detected_n, N_success)[2],
    rate_ci = format_pct_ci(detection_rate, ci_low, ci_high)
  ) |>
  dplyr::ungroup()

readr::write_csv(cell_summary, file.path(OUT, "data", "bernoulli_Q1_cell_summary.csv"))

gamma_summary <- cell_summary |>
  dplyr::group_by(gamma) |>
  dplyr::summarise(
    cells = dplyr::n(),
    total_success = sum(N_success),
    total_detected = sum(detected_n),
    mean_detection_rate = mean(detection_rate, na.rm = TRUE),
    pooled_detection_rate = ifelse(total_success > 0, total_detected / total_success, NA_real_),
    pooled_mc_se = mc_se(pooled_detection_rate, total_success),
    ci_low = wilson_ci(total_detected, total_success)[1],
    ci_high = wilson_ci(total_detected, total_success)[2],
    pooled_rate_ci = format_pct_ci(pooled_detection_rate, ci_low, ci_high),
    .groups = "drop"
  )

readr::write_csv(gamma_summary, file.path(OUT, "data", "bernoulli_Q1_summary_by_gamma.csv"))

size_summary <- cell_summary |>
  dplyr::select(m, n_i, gamma, detection_rate, mc_se, ci_low, ci_high, rate_ci) |>
  dplyr::arrange(m, n_i, gamma)

readr::write_csv(size_summary, file.path(OUT, "data", "bernoulli_Q1_size_summary_long.csv"))

best_stat_summary <- replicate_results |>
  dplyr::filter(gamma > 0, !is.na(best_stat)) |>
  dplyr::count(gamma, best_stat, sort = TRUE)

readr::write_csv(best_stat_summary, file.path(OUT, "data", "bernoulli_Q1_best_stat_summary.csv"))

# -----------------------------
# 9. LaTeX tables
# -----------------------------

latex_gamma <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{Bernoulli omitted-covariate power by omitted-effect magnitude. Detection is ",
         "based on the refit-calibrated minimum-$p$ omnibus statistic. Parentheses give ",
         "Wilson 95\\% Monte Carlo confidence intervals.}"),
  "\\label{tab:bernoulli-q1-power-gamma}",
  "\\begin{tabular}{rrrr}",
  "\\toprule",
  "$\\gamma$ & Total datasets & Detection rate & MC SE \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(gamma_summary))) {
  latex_gamma <- c(latex_gamma, sprintf(
    "%.2f & %d & %s & %.3f \\\\",
    gamma_summary$gamma[i],
    gamma_summary$total_success[i],
    gamma_summary$pooled_rate_ci[i],
    gamma_summary$pooled_mc_se[i]
  ))
}

latex_gamma <- c(latex_gamma, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(latex_gamma, file.path(OUT, "tables_latex", "table_bernoulli_Q1_power_by_gamma.tex"))

latex_cell <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{Bernoulli omitted-covariate power by design size. Detection rates include ",
         "Wilson 95\\% Monte Carlo confidence intervals.}"),
  "\\label{tab:bernoulli-q1-power-design}",
  "\\begin{tabular}{rrrrrr}",
  "\\toprule",
  "$m$ & $n_i$ & $\\gamma$ & Datasets & Detection rate & MC SE \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(cell_summary))) {
  latex_cell <- c(latex_cell, sprintf(
    "%d & %d & %.2f & %d & %s & %.3f \\\\",
    cell_summary$m[i],
    cell_summary$n_i[i],
    cell_summary$gamma[i],
    cell_summary$N_success[i],
    cell_summary$rate_ci[i],
    cell_summary$mc_se[i]
  ))
}

latex_cell <- c(latex_cell, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(latex_cell, file.path(OUT, "tables_latex", "table_bernoulli_Q1_power_by_design.tex"))

# -----------------------------
# 10. Figures
# -----------------------------

p1 <- ggplot2::ggplot(
  cell_summary,
  ggplot2::aes(x = gamma, y = detection_rate, group = factor(n_i), linetype = factor(n_i))
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.035) +
  ggplot2::facet_wrap(~ m, labeller = ggplot2::label_both) +
  ggplot2::theme_bw() +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(
    x = "Omitted covariate effect gamma",
    y = "Detection rate",
    linetype = "Cluster size n_i",
    title = "Bernoulli RR power curves with 95% Monte Carlo intervals"
  )

ggplot2::ggsave(
  file.path(OUT, "figures", "bernoulli_Q1_power_by_gamma_m_ni.png"),
  p1, width = 10, height = 6, dpi = 300
)

null_cells <- cell_summary |> dplyr::filter(gamma == 0)
upper_null <- suppressWarnings(max(null_cells$ci_high, na.rm = TRUE))
if (!is.finite(upper_null)) upper_null <- 0.20
upper_null <- max(0.20, upper_null)

p2 <- ggplot2::ggplot(
  null_cells,
  ggplot2::aes(x = factor(m), y = detection_rate, group = factor(n_i), linetype = factor(n_i))
) +
  ggplot2::geom_hline(yintercept = ALPHA, linetype = 2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.08) +
  ggplot2::theme_bw() +
  ggplot2::scale_y_continuous(limits = c(0, upper_null)) +
  ggplot2::labs(
    x = "Number of clusters m",
    y = "Null rejection rate",
    linetype = "Cluster size n_i",
    title = "Bernoulli RR null-level check with 95% Monte Carlo intervals"
  )

ggplot2::ggsave(
  file.path(OUT, "figures", "bernoulli_Q1_null_level_gamma0.png"),
  p2, width = 9, height = 5.5, dpi = 300
)

if (nrow(best_stat_summary) > 0) {
  p3 <- ggplot2::ggplot(
    best_stat_summary,
    ggplot2::aes(x = reorder(best_stat, n), y = n)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ gamma, scales = "free_y", labeller = ggplot2::label_both) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Most sensitive statistic",
      y = "Frequency",
      title = "Selected RR statistics under Bernoulli alternatives"
    )

  ggplot2::ggsave(
    file.path(OUT, "figures", "bernoulli_Q1_best_stat_frequency.png"),
    p3, width = 10, height = 6, dpi = 300
  )
}

# -----------------------------
# 11. Logs and manuscript note
# -----------------------------

run_end_time <- Sys.time()

method_note <- c(
  "RR_GLMM manuscript Bernoulli Q1-level power and null-level simulation",
  "",
  paste0("Run mode = ", RUN_MODE),
  paste0("N_ALT = ", N_ALT),
  paste0("N_NULL = ", N_NULL),
  paste0("B_CAL = ", B_CAL),
  paste0("N_QUAD = ", N_QUAD),
  paste0("m values = ", paste(M_VALUES, collapse = ", ")),
  paste0("n_i values = ", paste(N_VALUES, collapse = ", ")),
  paste0("gamma alternatives = ", paste(GAMMA_ALT, collapse = ", ")),
  paste0("gamma null = ", GAMMA_NULL),
  "",
  "Detection is based on p_combo_minP < 0.05.",
  "Monte Carlo SE and Wilson 95% CIs are reported for every design cell.",
  "",
  paste0("Run started: ", as.character(run_start_time)),
  paste0("Run ended: ", as.character(run_end_time)),
  paste0("Elapsed hours: ", round(as.numeric(difftime(run_end_time, run_start_time,
                                                      units = "hours")), 3)),
  paste0("Output folder: ", OUT)
)

writeLines(method_note, file.path(OUT, "logs", "README_manuscript_Bernoulli_Q1_power_level.txt"))

sink(file.path(OUT, "logs", "sessionInfo.txt"))
print(sessionInfo())
sink()

# -----------------------------
# 12. Completion message
# -----------------------------

log_msg("DONE.")
log_msg("Output folder: ", OUT)

cat("
============================================================
")
cat("DONE
")
cat("Run mode: ", RUN_MODE, "
", sep = "")
cat("Output folder:
", OUT, "
", sep = "")
cat("============================================================
")

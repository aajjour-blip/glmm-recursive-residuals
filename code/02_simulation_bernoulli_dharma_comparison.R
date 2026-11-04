# ============================================================
# RR_GLMM manuscript - R self-contained runner
# DHARMa head-to-head simulation for selected Bernoulli departures
#
# What this script does:
#   1. Loads required R packages.
#   2. Runs the proposed recursive-residual diagnostic vs DHARMa.
#   3. Saves replicate-level results, summaries, LaTeX tables, figures,
#      and logs into a structured output folder.
#
# Recommended R use:
#   source("code/02_simulation_bernoulli_dharma_comparison.R")
#
# Optional full-run overrides can be set with environment variables:
#   N_REP=200, B_CAL=200, DHARMA_NSIM=1000, M=100, N_I=8
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(stringsAsFactors = FALSE)
options(warn = 1)

# Parallel package compilation where possible.
if (requireNamespace("parallel", quietly = TRUE)) {
  ncpu <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  Sys.setenv(MAKEFLAGS = paste0("-j", ncpu))
}

is_R <- dir.exists(".") || nzchar(Sys.getenv("R_RELEASE_TAG"))

# Function definitions used in this script
source(file.path("code", "functions", "functions_02_simulation_bernoulli_dharma_comparison.R"), local = TRUE)

# System dependencies should be installed by the user before running this script.
required <- c("lme4", "statmod", "dplyr", "ggplot2", "readr", "tibble", "DHARMa")
load_required_packages(required)
suppressPackageStartupMessages({
  library(lme4)
  library(statmod)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})


FAST_TEST <- as_logical_env("FAST_TEST", FALSE)

set.seed(as_int_env("MASTER_SEED", 2026))

if (FAST_TEST) {
  # Smoke-test settings to confirm that the pipeline works before launching the full simulation.
  N_REP <- as_int_env("N_REP", 2)
  B_CAL <- as_int_env("B_CAL", 50)  # >= 40 required by RR calibration
  DHARMA_NSIM <- as_int_env("DHARMA_NSIM", 100)
  M <- as_int_env("M", 30)
  N_I <- as_int_env("N_I", 5)
} else {
  # Manuscript-strength settings. These can be overridden using environment variables.
  N_REP <- as_int_env("N_REP", 200)
  B_CAL <- as_int_env("B_CAL", 200)
  DHARMA_NSIM <- as_int_env("DHARMA_NSIM", 1000)
  M <- as_int_env("M", 100)
  N_I <- as_int_env("N_I", 8)
}

ALPHA <- as_num_env("ALPHA", 0.05)
N_QUAD <- as_int_env("N_QUAD", 25)

BETA0 <- as_num_env("BETA0", -1.00)
BETA1 <- as_num_env("BETA1",  0.80)
SIGMA_B <- as_num_env("SIGMA_B", 0.80)

SCENARIOS <- c("null", "omitted_x2_gamma05", "omitted_x2_gamma10", "serial_dependence_rho10")
RESUME <- as_logical_env("RESUME", TRUE)

OUT_ROOT <- Sys.getenv("OUT_ROOT", unset = file.path(getwd(), "results"))
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
OUT_NAME <- Sys.getenv("OUT_NAME", unset = "bernoulli_dharma_comparison_results")
OUT <- file.path(OUT_ROOT, OUT_NAME)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tables_latex"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)

run_start_time <- Sys.time()

cat("============================================================\n")
cat("RR_GLMM manuscript DHARMa head-to-head binary simulation\n")
cat("R mode:", is_R, "\n")
cat("FAST_TEST:", FAST_TEST, "\n")
cat("N_REP:", N_REP, "| B_CAL:", B_CAL, "| DHARMA_NSIM:", DHARMA_NSIM, "\n")
cat("M:", M, "| N_I:", N_I, "| N_QUAD:", N_QUAD, "\n")
cat("Output folder:", OUT, "\n")
cat("============================================================\n\n")

writeLines(c(
  "RR_GLMM manuscript DHARMa head-to-head binary comparison - R runner",
  paste0("Output folder: ", OUT),
  "",
  "To run a smoke test first:",
  "  Sys.setenv(FAST_TEST = TRUE)",
  "  source(\"code/02_simulation_bernoulli_dharma_comparison.R\")"
), file.path(OUT, "logs", "run_summary.txt"))

# ============================================================
# Reuse utility functions from the main manuscript script
# ============================================================


# ============================================================
# DGP
# ============================================================


# ============================================================
# Run grid with resume
# ============================================================

design_grid <- expand.grid(
  scenario = SCENARIOS,
  rep_id = seq_len(N_REP),
  KEEP.OUT.ATTRS = FALSE
) |> as_tibble() |> arrange(scenario, rep_id)

write_csv(design_grid, file.path(OUT, "data", "dharma_head_to_head_design_grid.csv"))

partial_file <- file.path(OUT, "data", "dharma_head_to_head_replicate_results_partial.csv")

completed <- tibble()
if (RESUME && file.exists(partial_file)) {
  completed <- tryCatch(read_csv(partial_file, show_col_types = FALSE),
                        error = function(e) tibble())
  if (nrow(completed) > 0) {
    completed <- completed |> distinct(scenario, rep_id, .keep_all = TRUE)
  }
}

if (nrow(completed) > 0) {
  remaining <- design_grid |> anti_join(completed |> select(scenario, rep_id),
                                        by = c("scenario", "rep_id"))
} else {
  remaining <- design_grid
}

cat("Total rows:", nrow(design_grid), "\n")
cat("Completed rows:", nrow(completed), "\n")
cat("Remaining rows:", nrow(remaining), "\n")

new_results <- vector("list", nrow(remaining))

if (nrow(remaining) > 0) {
  for (ii in seq_len(nrow(remaining))) {
    row <- remaining[ii, ]
    elapsed <- difftime(Sys.time(), run_start_time, units = "mins")

    cat(
      "Row", ii, "of", nrow(remaining),
      "| scenario =", row$scenario,
      "| rep =", row$rep_id,
      "| elapsed min =", round(as.numeric(elapsed), 2), "\n"
    )

    res <- run_one_comparison(row$scenario, row$rep_id)
    new_results[[ii]] <- res

    if (!file.exists(partial_file) || file.info(partial_file)$size == 0) {
      write_csv(res, partial_file)
    } else {
      write_csv(res, partial_file, append = TRUE)
    }

    if (ii %% 5 == 0) gc(verbose = FALSE)
  }
}

new_df <- bind_rows(new_results)

replicate_results <- bind_rows(completed, new_df) |>
  distinct(scenario, rep_id, .keep_all = TRUE) |>
  arrange(scenario, rep_id)

write_csv(replicate_results, file.path(OUT, "data", "dharma_head_to_head_replicate_results.csv"))

# ============================================================
# Summaries with MC SE and Wilson CI
# ============================================================


s_rr <- summarise_detector(replicate_results, "rr_detect")
s_raw <- summarise_detector(replicate_results, "dharma_any_raw_detect")
s_bonf <- summarise_detector(replicate_results, "dharma_any_bonf_detect")

summary <- bind_rows(s_rr, s_raw, s_bonf)

write_csv(summary, file.path(OUT, "data", "dharma_head_to_head_detection_summary.csv"))

best_rr <- replicate_results |>
  filter(!is.na(rr_best_stat)) |>
  count(scenario, rr_best_stat, sort = TRUE)

write_csv(best_rr, file.path(OUT, "data", "dharma_head_to_head_rr_best_stat_summary.csv"))

# LaTeX table
latex <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{Head-to-head simulation comparison between recursive residuals and DHARMa for ",
         "selected Bernoulli departures. Detection rates include Wilson 95\\% Monte Carlo ",
         "confidence intervals.}"),
  "\\label{tab:dharma-head-to-head-bernoulli}",
  "\\begin{tabular}{llrrr}",
  "\\toprule",
  "Scenario & Detector & Datasets & Detection rate & MC SE \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(summary))) {
  latex <- c(latex, sprintf(
    "%s & %s & %d & %.1f\\%% (%.1f, %.1f) & %.3f \\\\",
    summary$scenario[i],
    summary$detector[i],
    summary$N_success[i],
    100 * summary$detection_rate[i],
    100 * summary$ci_low[i],
    100 * summary$ci_high[i],
    summary$mc_se[i]
  ))
}

latex <- c(latex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(latex, file.path(OUT, "tables_latex", "table_DHARMa_head_to_head_bernoulli.tex"))

# Figure
p <- ggplot(summary, aes(x = scenario, y = detection_rate, fill = detector)) +
  geom_col(position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.75), width = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Scenario",
    y = "Detection rate",
    fill = "Detector",
    title = "RR vs DHARMa detection in Bernoulli GLMM simulations"
  )

ggsave(file.path(OUT, "figures", "RR_vs_DHARMa_bernoulli_detection.png"),
       p, width = 9, height = 5.5, dpi = 300)

# Logs
run_end_time <- Sys.time()

writeLines(c(
  "RR_GLMM manuscript DHARMa head-to-head binary comparison",
  paste0("N_REP = ", N_REP),
  paste0("B_CAL = ", B_CAL),
  paste0("DHARMA_NSIM = ", DHARMA_NSIM),
  paste0("m = ", M),
  paste0("n_i = ", N_I),
  paste0("Run started: ", as.character(run_start_time)),
  paste0("Run ended: ", as.character(run_end_time)),
  paste0("Elapsed hours: ", round(as.numeric(difftime(run_end_time, run_start_time,
                                                      units = 'hours')), 3))
), file.path(OUT, "logs", "README_DHARMa_head_to_head_binary.txt"))

sink(file.path(OUT, "logs", "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("
DONE. Results were written to the output folder.
")
cat("Output folder:
", OUT, "
")

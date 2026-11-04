# CBPP empirical example summary and verification.
#
# Purpose:
#   Recreate the manuscript CBPP table from the archived B = 1000 recursive-
#   residual and DHARMa/AIC output files.
#
# Arguments:
#   None. Run from the root of the code and data supplement.
#
# Returns:
#   Writes refreshed Table 8 to tables/ and a verification CSV to
#   results/empirical_cbpp/.

root <- getwd()
out_dir <- file.path(root, "results", "empirical_cbpp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- file.path(root, "data", "empirical")
tab_dir <- file.path(root, "tables")
rr_file <- file.path(data_dir, "cbpp_RR_B1000_component_calibration.csv")
dharma_file <- file.path(data_dir, "cbpp_DHARMa_AIC_B1000.csv")
run_file <- file.path(data_dir, "cbpp_run_summary.csv")

for (f in c(rr_file, dharma_file, run_file)) if (!file.exists(f)) stop("Missing: ", f)

rr <- read.csv(rr_file, stringsAsFactors = FALSE)
dh <- read.csv(dharma_file, stringsAsFactors = FALSE)
run <- read.csv(run_file, stringsAsFactors = FALSE)

required_rr <- c("p_combo_minP", "p_min_raw", "best_stat", "detected_combo")
missing_rr <- setdiff(required_rr, names(rr))
if (length(missing_rr)) stop("Missing columns in RR file: ", paste(missing_rr, collapse = ", "))
required_dh <- c("DHARMa_uniformity_p", "DHARMa_dispersion_p", "AIC_baseline", "AIC_OLRE")
missing_dh <- setdiff(required_dh, names(dh))
if (length(missing_dh)) {
  stop("Missing columns in DHARMa/AIC file: ", paste(missing_dh, collapse = ", "))
}
required_run <- c("B_ENV", "N_QUAD", "cbpp_refit_failures", "cbpp_singular_refits")
missing_run <- setdiff(required_run, names(run))
if (length(missing_run)) {
  stop("Missing columns in run summary: ", paste(missing_run, collapse = ", "))
}

fmt_stat <- function(z) {
  map <- c(
    Tinf = "Tinf", T2 = "T2", Tmax_r = "Tmax_r", Tmax_C = "Tmax_C",
    Tdisp = "Tdisp", Tlag = "Tlag"
  )
  ifelse(z %in% names(map), map[z], z)
}

p_combo <- rr$p_combo_minP[1]
p_min <- rr$p_min_raw[1]
best <- fmt_stat(rr$best_stat[1])
reject <- ifelse(isTRUE(rr$detected_combo[1]), "yes", "no")

lines <- c(
  "\\begin{table}[htbp]", "\\centering", "\\small",
  "\\caption{CBPP clustered-binomial empirical diagnostic summary.}",
  "\\label{tab:empirical-cbpp}", "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  paste0(
    "Model & RR $p_{\\mathrm{combo}}$ & RR $p_{\\min}$ & Best RR & ",
    "RR reject & \\texttt{DHARMa} unif. & \\texttt{DHARMa} disp. & AIC \\\\"
  ),
  "\\midrule",
  sprintf(
    "Binomial RI, period + herd & %.6f & %.6f & %s & %s & %.3f & %.3f & %.1f \\\\",
    p_combo, p_min, best, reject, dh$DHARMa_uniformity_p[1], dh$DHARMa_dispersion_p[1],
    dh$AIC_baseline[1]
  ),
  sprintf("Binomial RI + OLRE & --$^{\\dag}$ & -- & -- & -- & -- & -- & %.1f \\\\", dh$AIC_OLRE[1]),
  "\\bottomrule", "\\end{tabular}}%",
  paste0(
    "\\TableNote{$^{\\dag}$The RR diagnostic was applied to the simpler binomial ",
    "random-intercept model only; the OLRE model is listed for AIC comparison. ",
    "The RR calibration used $B=", run$B_ENV[1], "$ refit replicates and ",
    "$N_{\\mathrm{quad}}=", run$N_QUAD[1], "$. The CBPP refit calibration had ",
    run$cbpp_singular_refits[1], " singular simulated refits and ",
    run$cbpp_refit_failures[1], " refit failures; singular simulated refits are part ",
    "of the fitted-model calibration distribution rather than a failure of the observed model.}"
  ),
  "\\end{table}"
)
writeLines(lines, file.path(tab_dir, "table_empirical_cbpp.tex"))

summary <- data.frame(
  p_combo = p_combo,
  p_min = p_min,
  best_stat = best,
  reject = reject,
  dharma_uniformity_p = dh$DHARMa_uniformity_p[1],
  dharma_dispersion_p = dh$DHARMa_dispersion_p[1],
  AIC_baseline = dh$AIC_baseline[1],
  AIC_OLRE = dh$AIC_OLRE[1],
  B = run$B_ENV[1],
  Nquad = run$N_QUAD[1],
  singular_refits = run$cbpp_singular_refits[1],
  refit_failures = run$cbpp_refit_failures[1],
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "cbpp_verification.csv"), row.names = FALSE)
cat("CBPP empirical table refreshed and verified.\n")
print(summary)

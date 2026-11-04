# Simulation catalogue summary and verification.
#
# Purpose:
#   Recreate the manuscript and Supporting Information summaries for the
#   non-Gaussian simulation catalogue from the archived B = 1000 output table.
#   The heavy scenario-level simulations are archived in data/simulation_catalogue/;
#   this script formats and verifies the reported tables from those outputs.
#
# Arguments:
#   None. Run from the root of the code and data supplement.
#
# Returns:
#   Writes refreshed catalogue tables to tables/ and verification summaries to
#   results/simulation_catalogue/.

root <- getwd()
out_dir <- file.path(root, "results", "simulation_catalogue")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- file.path(root, "data", "simulation_catalogue")
tab_dir <- file.path(root, "tables")
primary_file <- file.path(data_dir, "simulation_catalogue_primary_B1000.csv")
family_file <- file.path(data_dir, "family_detection_rates.csv")

if (!file.exists(primary_file)) stop("Missing: ", primary_file)
if (!file.exists(family_file)) stop("Missing: ", family_file)

x <- read.csv(primary_file, stringsAsFactors = FALSE)
fam <- read.csv(family_file, stringsAsFactors = FALSE)

required <- c(
  "scenario", "family", "q_label", "p_Tinf", "p_T2", "p_Tmax_r",
  "p_Tmax_C", "p_Tdisp", "p_Tjump", "p_Tlag", "p_min_component",
  "p_combo_omnibus", "best_component_statistic", "detected_combo_omnibus"
)
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing columns in catalogue file: ", paste(missing, collapse = ", "))

# Python-produced archive encodes booleans as text "True"/"False".
# R read.csv() does not consistently interpret these as logical values.
# Convert strictly, with no coercion of missing or unexpected values.
raw_detected <- trimws(as.character(x$detected_combo_omnibus))
if (anyNA(raw_detected) || !all(raw_detected %in% c("True", "False", "TRUE", "FALSE"))) {
  stop("Invalid detected_combo_omnibus values; expected True/False only")
}
x$detected_combo_omnibus <- raw_detected %in% c("True", "TRUE")

# Confirm that archived binary decisions agree with the stated 5% rule.
if (anyNA(x$p_combo_omnibus) ||
    any(x$detected_combo_omnibus != (x$p_combo_omnibus <= 0.05))) {
  stop("Archived omnibus decisions disagree with p_combo_omnibus <= 0.05")
}

# Verify family summary against underlying scenario records before formatting.
for (j in seq_len(nrow(fam))) {
  these <- x$family == fam$family[j]
  if (anyNA(these) || sum(these) != fam$scenarios[j] ||
      sum(x$detected_combo_omnibus[these]) != fam$detected[j] ||
      abs(mean(x$detected_combo_omnibus[these]) - fam$rate[j]) > 1e-9) {
    stop("Family detection summary disagrees with primary CSV: ", fam$family[j])
  }
}

escape_tex <- function(z) {
  z <- gsub("_", "\\_", z, fixed = TRUE)
  z
}

family_short <- function(z) {
  out <- ifelse(z == "binomial", "Bin.", ifelse(z == "poisson", "Pois.", "Gam."))
  out
}

stat_tex <- function(z) {
  map <- c(
    Tinf = "$T_\\infty$", T2 = "$T_2$", Tmax_r = "$T_{\\max r}$",
    Tmax_C = "$T_{\\max C}$", Tdisp = "$T_{\\mathrm{disp}}$",
    Tjump = "$T_{\\mathrm{jump}}$", Tlag = "$T_{\\mathrm{lag}}$"
  )
  unname(ifelse(z %in% names(map), map[z], z))
}

fmt_p <- function(z) sprintf("%.3f", z)

# Table 3: detection rates by family.
fam_order <- c("binomial", "poisson", "gamma")
fam <- fam[match(fam_order, fam$family), ]
fam_lab <- c(binomial = "Bernoulli", poisson = "Poisson", gamma = "Gamma")
lines <- c(
  "\\begin{table}[htbp]",
  "\\small",
  "\\caption{Detection rates by response family using the refit-calibrated omnibus",
  "statistic at the $5\\%$ level with $B = 1000$ replicates.",
  "A scenario is counted as detected if $p_{\\mathrm{combo}} \\leq 0.05$ in that run.",
  "Correctly specified models from each family were not rejected.}",
  "\\label{tab:target-family-detection}",
  "\\begin{tabular*}{\\textwidth}{@{}l@{\\extracolsep{\\fill}}rrr@{}}",
  "\\toprule",
  "Family & Scenarios Detected & Total Scenarios & Detection Rate \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(fam))) {
  lines <- c(lines, sprintf(
    "%s & %d & %d & %.1f\\%% \\\\", fam_lab[[fam$family[i]]], fam$detected[i],
    fam$scenarios[i], 100 * fam$rate[i]
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular*}", "\\end{table}")
writeLines(lines, file.path(tab_dir, "table_detection_by_family.tex"))

# Supporting Information: full catalogue table.
lines <- c(
  "\\begingroup",
  "\\tiny",
  "\\setlength{\\tabcolsep}{1.8pt}",
  "\\renewcommand{\\arraystretch}{0.92}",
  "\\begin{longtable}{p{0.30\\linewidth}p{0.05\\linewidth}p{0.20\\linewidth}rrrrrl}",
  paste0(
    "\\caption{Full calibrated recursive-residual diagnostic results for all ",
    "non-Gaussian simulation scenarios using $B=1000$ calibration replicates. ",
    "Columns $p_\\infty$ and $p_2$ report p-values for $T_\\infty$ and $T_2$. ",
    "The formal global decision uses $p_{\\mathrm{combo}}$; the raw $p_{\\min}$ ",
    "is retained as a screening measure.}\\label{tab:simulation-catalogue-full}\\\\"
  ),
  "\\toprule",
  paste0(
    "Scenario & Fam. & Index & $p_\\infty$ & $p_2$ & $p_{\\min}$ & ",
    "$p_{\\mathrm{combo}}$ & Best & Det. \\\\"
  ),
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  paste0(
    "Scenario & Fam. & Index & $p_\\infty$ & $p_2$ & $p_{\\min}$ & ",
    "$p_{\\mathrm{combo}}$ & Best & Det. \\\\"
  ),
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot"
)
for (i in seq_len(nrow(x))) {
  lines <- c(lines, sprintf(
    "%s & %s & %s & %.3f & %.3f & %.3f & %.3f & %s & %s \\\\",
    escape_tex(x$scenario[i]), family_short(x$family[i]), escape_tex(x$q_label[i]),
    x$p_Tinf[i], x$p_T2[i], x$p_min_component[i], x$p_combo_omnibus[i],
    stat_tex(x$best_component_statistic[i]), ifelse(x$detected_combo_omnibus[i], "Y", "N")
  ))
}
lines <- c(lines, "\\end{longtable}", "\\endgroup")
writeLines(lines, file.path(tab_dir, "table_simulation_catalogue_full_results.tex"))

# Detected and borderline tables.
det <- x[x$detected_combo_omnibus, ]
lines <- c(
  "\\begin{table}[H]", "\\centering", "\\scriptsize",
  paste0(
    "\\caption{Scenarios detected by the calibrated omnibus recursive-residual ",
    "statistic in the final $B=1000$ run.}"
  ),
  "\\label{tab:detected-scenarios}", "\\begin{tabular}{lllr}", "\\toprule",
  "Scenario & Family & Best statistic & $p_{\\mathrm{combo}}$ \\\\", "\\midrule"
)
for (i in seq_len(nrow(det))) {
  lines <- c(lines, sprintf(
    "%s & %s & %s & %.3f \\\\", escape_tex(det$scenario[i]), family_short(det$family[i]),
    stat_tex(det$best_component_statistic[i]), det$p_combo_omnibus[i]
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(tab_dir, "table_detected_scenarios.tex"))

border <- x[!x$detected_combo_omnibus & x$p_min_component <= 0.05, ]
border <- border[order(border$p_combo_omnibus), ]
lines <- c(
  "\\begin{table}[H]", "\\centering", "\\scriptsize",
  "\\caption{Borderline and screening-positive scenarios in the final $B=1000$ run.}",
  "\\label{tab:borderline-scenarios}", "\\begin{tabular}{llrrl}", "\\toprule",
  "Scenario & Family & $p_{\\min}$ & $p_{\\mathrm{combo}}$ & Interpretation \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(border))) {
  lines <- c(lines, sprintf(
    "%s & %s & %.3f & %.3f & screening only \\\\", escape_tex(border$scenario[i]),
    family_short(border$family[i]), border$p_min_component[i], border$p_combo_omnibus[i]
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(tab_dir, "table_borderline_scenarios.tex"))

summary <- data.frame(
  n_scenarios = nrow(x),
  n_detected = sum(x$detected_combo_omnibus),
  n_borderline_screening = nrow(border),
  n_best_Tjump = sum(x$best_component_statistic == "Tjump"),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "simulation_catalogue_verification.csv"), row.names = FALSE)
cat("Simulation catalogue summary refreshed and verified.\n")
print(summary)

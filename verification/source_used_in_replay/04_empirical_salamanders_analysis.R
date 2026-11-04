# ============================================================
# Salamanders empirical RR diagnostics, B = 1000
# Poisson RI, NB2 RI, and ZINB RI comparison
#
# Run from the root of the code/data supplement with:
#   source("code/04_empirical_salamanders_analysis.R")
#
# Outputs are written under ./results/RR_Salamanders_RR_B1000_NB2_ZINB_submission_results.
# ============================================================

required <- c("glmmTMB", "DHARMa", "statmod", "dplyr", "ggplot2", "readr", "tibble")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install missing packages first: ", paste(missing, collapse = ", "))
}

library(glmmTMB)
library(DHARMa)
library(statmod)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(2026)

B_ENV <- 1000
N_QUAD <- 50
ALPHA <- 0.05

# Switch models on/off here if runtime becomes too long.
RUN_POISSON <- TRUE
RUN_NB2     <- TRUE
RUN_ZINB    <- TRUE

# Function definitions used in this script
source(file.path("code", "functions", "functions_04_empirical_salamanders_analysis.R"), local = TRUE)

OUT <- file.path(get_output_root(), "RR_Salamanders_RR_B1000_NB2_ZINB_submission_results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "logs"), showWarnings = FALSE, recursive = TRUE)

cat("\nOutputs will be written to:\n", OUT, "\n\n")

# ============================================================
# General utility functions
# ============================================================


# ============================================================
# Salamanders model-specific helper functions
# ============================================================


# ============================================================
# Recursive residuals for Salamanders random-intercept count GLMMs
#
# Supports:
#   family_type = "poisson", zero_inflated = FALSE
#   family_type = "nbinom2", zero_inflated = FALSE
#   family_type = "nbinom2", zero_inflated = TRUE
#
# Baseline fixed-effects formula:
#   count ~ spp * mined + (1 | site)
#
# ZI formula for ZINB:
#   ziformula = ~ spp * mined
# ============================================================


# ============================================================
# Fit Salamanders models
# ============================================================

data(Salamanders, package = "glmmTMB")

Salamanders <- Salamanders |>
  mutate(
    site = factor(site),
    sample_num = suppressWarnings(as.numeric(as.character(sample))),
    mined = factor(mined),
    spp = factor(spp)
  )

if (any(is.na(Salamanders$sample_num))) {
  Salamanders$sample_num <- as.numeric(Salamanders$sample)
}

fit_sal_pois <- glmmTMB(count ~ spp * mined + (1 | site),
                        data = Salamanders,
                        family = poisson)

fit_sal_nb2 <- glmmTMB(count ~ spp * mined + (1 | site),
                       data = Salamanders,
                       family = nbinom2)

fit_sal_zip <- glmmTMB(count ~ spp * mined + (1 | site),
                       ziformula = ~ spp * mined,
                       data = Salamanders,
                       family = poisson)

fit_sal_zinb <- glmmTMB(count ~ spp * mined + (1 | site),
                        ziformula = ~ spp * mined,
                        data = Salamanders,
                        family = nbinom2)

sal_aic <- AIC(fit_sal_pois, fit_sal_nb2, fit_sal_zip, fit_sal_zinb) |>
  as.data.frame() |>
  tibble::rownames_to_column("model")

readr::write_csv(sal_aic, file.path(OUT, "data", "salamanders_model_AIC.csv"))

# ============================================================
# Function to run RR calibration for one Salamanders model
# ============================================================


# ============================================================
# Run selected RR models
# ============================================================

summary_list <- list()

if (RUN_POISSON) {
  summary_list[["Poisson_RI"]] <- run_rr_model(
    model_name = "Poisson_RI",
    fit = fit_sal_pois,
    family_type = "poisson",
    zero_inflated = FALSE
  )
}

if (RUN_NB2) {
  summary_list[["NB2_RI"]] <- run_rr_model(
    model_name = "NB2_RI",
    fit = fit_sal_nb2,
    family_type = "nbinom2",
    zero_inflated = FALSE
  )
}

if (RUN_ZINB) {
  summary_list[["ZINB_RI"]] <- run_rr_model(
    model_name = "ZINB_RI",
    fit = fit_sal_zinb,
    family_type = "nbinom2",
    zero_inflated = TRUE
  )
}

empirical_rr_summary <- dplyr::bind_rows(summary_list)
readr::write_csv(empirical_rr_summary, file.path(OUT, "data",
                                                 "salamanders_RR_all_models_summary_B1000.csv"))

# DHARMa is deliberately run after RR calibration to avoid comparator-induced
# cross-version effects on the RR reference simulation.
# ============================================================
# DHARMa diagnostics for modern comparators
# ============================================================

cat("Running DHARMa diagnostics...\n")

sim_sal_pois <- DHARMa::simulateResiduals(fit_sal_pois, n = B_ENV, seed = 3026)
sim_sal_nb2  <- DHARMa::simulateResiduals(fit_sal_nb2,  n = B_ENV, seed = 3027)
sim_sal_zip  <- DHARMa::simulateResiduals(fit_sal_zip,  n = B_ENV, seed = 3028)
sim_sal_zinb <- DHARMa::simulateResiduals(fit_sal_zinb, n = B_ENV, seed = 3029)

png(file.path(OUT, "figures", "salamanders_DHARMa_poisson.png"), 1200, 900)
plot(sim_sal_pois)
dev.off()

png(file.path(OUT, "figures", "salamanders_DHARMa_nb2.png"), 1200, 900)
plot(sim_sal_nb2)
dev.off()

png(file.path(OUT, "figures", "salamanders_DHARMa_zinb.png"), 1200, 900)
plot(sim_sal_zinb)
dev.off()

sal_dharma <- data.frame(
  dataset = "Salamanders",
  model = c("Poisson RI", "NB2 RI", "ZIP RI", "ZINB RI"),
  DHARMa_uniformity_p = c(testUniformity(sim_sal_pois)$p.value,
                          testUniformity(sim_sal_nb2)$p.value,
                          testUniformity(sim_sal_zip)$p.value,
                          testUniformity(sim_sal_zinb)$p.value),
  DHARMa_dispersion_p = c(testDispersion(sim_sal_pois)$p.value,
                          testDispersion(sim_sal_nb2)$p.value,
                          testDispersion(sim_sal_zip)$p.value,
                          testDispersion(sim_sal_zinb)$p.value),
  DHARMa_zeroInflation_p = c(testZeroInflation(sim_sal_pois)$p.value,
                             testZeroInflation(sim_sal_nb2)$p.value,
                             testZeroInflation(sim_sal_zip)$p.value,
                             testZeroInflation(sim_sal_zinb)$p.value),
  DHARMa_outlier_p = c(testOutliers(sim_sal_pois, type = "bootstrap")$p.value,
                       testOutliers(sim_sal_nb2, type = "bootstrap")$p.value,
                       testOutliers(sim_sal_zip, type = "bootstrap")$p.value,
                       testOutliers(sim_sal_zinb, type = "bootstrap")$p.value)
)

readr::write_csv(sal_dharma, file.path(OUT, "data", "salamanders_DHARMa_tests_B1000.csv"))


# ============================================================
# Final logs
# ============================================================

sink(file.path(OUT, "logs", "sessionInfo_salamanders_RR_B1000_NB2_ZINB.txt"))
print(sessionInfo())
sink()

cat("\nDONE.\n")
cat("All outputs are in:\n")
cat(OUT, "\n")
cat("\nPlease zip this folder and save the output:\n")
cat("RR_Salamanders_RR_B1000_NB2_ZINB_submission_results\n")

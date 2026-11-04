# ============================================================
# RR-guided GLMM model refinement, B = 1000
# Run from the root of the code/data supplement with:
#   source("code/03_empirical_epilepsy_refinement.R")
# Outputs are written under ./results/RR_guided_model_refinement_B1000_submission_results
# ============================================================

required <- c("glmmTMB", "DHARMa", "statmod", "MASS", "dplyr", "ggplot2", "readr", "tibble")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) stop("Please install missing packages first: ", paste(missing,
    collapse = ", "))

library(glmmTMB)
library(DHARMa)
library(statmod)
library(MASS)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(2026)
B_ENV <- 1000
N_QUAD <- 40
ALPHA <- 0.05

RUN_SALAMANDERS_NB2_ENV <- TRUE
RUN_OWLS_ZINB_ARRIVAL  <- TRUE
RUN_EPIL_NB2_REFINED   <- TRUE

# Function definitions used in this script
source(file.path("code", "functions", "functions_03_empirical_epilepsy_refinement.R"), local = TRUE)

OUT <- file.path(get_output_root(), "RR_guided_model_refinement_B1000_submission_results")
dir.create(file.path(OUT, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
cat("\nOutputs will be written to:\n", OUT, "\n\n")

previous_reference <- data.frame(
  example = c("Salamanders", "Owls", "epil"),
  previous_model = c("NB2_RI", "ZINB_RI", "NB2_RI"),
  previous_RR_p_combo = c(0.416583, 0.000999, 0.000999),
  previous_best_stat = c("Tlag", "Tinf_ArrivalTime_num", "T2_base"),
  previous_interpretation = c(
    "NB2 already removed the main Poisson RR lack-of-fit signal.",
    "ZINB still rejected; RR suggested residual arrival-time/process structure.",
    "NB2 still rejected; RR suggested residual baseline/period structure."
  )
)
write_csv(previous_reference, file.path(OUT, "data", "previous_RR_reference_for_refinement.csv"))

# ---------------- Utility functions ----------------


# ---------------- Generic RR for count RI GLMMs ----------------


all_rr_summaries <- list(); all_dharma <- list(); all_aic <- list()

# ---------------- 1. Salamanders: NB2 with environmental/time covariates ----------------
if (RUN_SALAMANDERS_NB2_ENV) {
  cat("\nLoading Salamanders refinement...\n")
  data(Salamanders, package = "glmmTMB")
  sal <- Salamanders |>
    mutate(site = factor(site), sample_num = suppressWarnings(as.numeric(as.character(sample))),
           mined = factor(mined), spp = factor(spp), cover_s = as.numeric(scale(cover)),
           DOP_s = as.numeric(scale(DOP)), Wtemp_s = as.numeric(scale(Wtemp)),
           DOY_s = as.numeric(scale(DOY)), sample_s = as.numeric(scale(sample_num)))
  if (any(is.na(sal$sample_num))) { sal$sample_num <- as.numeric(sal$sample)
    sal$sample_s <- as.numeric(scale(sal$sample_num)) }

  fit_base <- glmmTMB(count ~ spp * mined + (1 | site), data = sal, family = nbinom2)
  fit_env <- glmmTMB(count ~ spp * mined + cover_s + DOP_s + Wtemp_s + DOY_s + sample_s +
                     (1 | site), data = sal, family = nbinom2)
  aic <- AIC(fit_base, fit_env) |> as.data.frame() |> tibble::rownames_to_column("model")
  write_csv(aic, file.path(OUT, "data",
                           "Salamanders_refinement_AIC.csv"))
                           all_aic[["Salamanders"]] <- data.frame(example = "Salamanders", aic)
  fixed <- ~ spp * mined + cover_s + DOP_s + Wtemp_s + DOY_s + sample_s
  refit <- function(d) glmmTMB(count ~ spp * mined + cover_s + DOP_s + Wtemp_s + DOY_s +
                               sample_s + (1 | site), data = d, family = nbinom2)
  all_rr_summaries[["Salamanders_NB2_env"]] <- run_count_rr_model("Salamanders", "NB2_env",
      fit_env, sal, fixed, refit, "site", "sample_num", "count", "nbinom2", FALSE, NULL, NULL,
          c("sample_num", "pred_mean", "cover", "DOP", "Wtemp", "DOY"))
}

# ---------------- 2. Owls: ZINB with arrival-time terms ----------------
if (RUN_OWLS_ZINB_ARRIVAL) {
  cat("\nLoading Owls refinement...\n")
  data(Owls, package = "glmmTMB")
  owls <- Owls |>
    mutate(Nest = factor(Nest), FoodTreatment = factor(FoodTreatment),
           SexParent = factor(SexParent),
           ArrivalTime_num = suppressWarnings(as.numeric(ArrivalTime)),
               ArrivalTime_s = as.numeric(scale(ArrivalTime_num)),
           logBroodSize = log(BroodSize))
  fit_base <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + offset(logBroodSize) +
                      (1 | Nest),
                      ziformula = ~ FoodTreatment * SexParent, data = owls, family = nbinom2)
  fit_arr <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent + poly(ArrivalTime_s,
      2) + offset(logBroodSize) + (1 | Nest),
                     ziformula = ~ FoodTreatment * SexParent + poly(ArrivalTime_s, 2),
                         data = owls, family = nbinom2)
  aic <- AIC(fit_base, fit_arr) |> as.data.frame() |> tibble::rownames_to_column("model")
  write_csv(aic, file.path(OUT, "data",
                           "Owls_refinement_AIC.csv"))
                           all_aic[["Owls"]] <- data.frame(example = "Owls", aic)
  fixed <- ~ FoodTreatment * SexParent + poly(ArrivalTime_s, 2)
  zi <- ~ FoodTreatment * SexParent + poly(ArrivalTime_s, 2)
  rr_offset_function <- function(d) d$logBroodSize
  refit <- function(d) glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                               poly(ArrivalTime_s, 2) + offset(logBroodSize) + (1 | Nest),
                              ziformula = ~ FoodTreatment * SexParent + poly(ArrivalTime_s, 2),
                               data = d, family = nbinom2)
  all_rr_summaries[["Owls_ZINB_arrival"]] <- run_count_rr_model("Owls", "ZINB_arrival",
      fit_arr, owls, fixed, refit, "Nest", "ArrivalTime_num", "SiblingNegotiation", "nbinom2",
          TRUE, zi, rr_offset_function, c("ArrivalTime_num", "pred_mean", "FoodTreatment", "SexParent",
                              "BroodSize"))
}

# ---------------- 3. epil: NB2 with nonlinear baseline + treatment-period ----------------
if (RUN_EPIL_NB2_REFINED) {
  cat("\nLoading epil refinement...\n")
  data(epil, package = "MASS")
  ep <- epil
  if (!"subject" %in% names(ep)) { if ("id" %in% names(ep)) ep$subject <- ep$id
    if ("ID" %in% names(ep)) ep$subject <- ep$ID }
  ep <- ep |>
    mutate(subject = factor(subject),
           period_num = suppressWarnings(as.numeric(as.character(period))),
           period_fac = factor(period), trt = factor(trt), log_base = log(base + 1),
           log_base_s = as.numeric(scale(log_base)), age_s = as.numeric(scale(age)))
  if (any(is.na(ep$period_num))) ep$period_num <- as.numeric(ep$period)

  fit_base <- glmmTMB(y ~ trt + base + age + period_fac + (1 | subject), data = ep,
                      family = nbinom2)
  fit_ref <- glmmTMB(y ~ trt * period_fac + poly(log_base_s, 2) + age_s + (1 | subject),
                     data = ep, family = nbinom2)
  aic <- AIC(fit_base, fit_ref) |> as.data.frame() |> tibble::rownames_to_column("model")
  write_csv(aic, file.path(OUT, "data",
                           "epil_refinement_AIC.csv"))
                           all_aic[["epil"]] <- data.frame(example = "epil", aic)
  fixed <- ~ trt * period_fac + poly(log_base_s, 2) + age_s
  refit <- function(d) glmmTMB(y ~ trt * period_fac + poly(log_base_s,
                                                           2) + age_s + (1 | subject),
                                                               data = d, family = nbinom2)
  all_rr_summaries[["epil_NB2_refined"]] <- run_count_rr_model("epil", "NB2_refined", fit_ref,
      ep, fixed, refit, "subject", "period_num", "y", "nbinom2", FALSE, NULL, NULL,
          c("period_num", "pred_mean", "base", "log_base", "age", "trt"))
}

# Run comparator diagnostics only after RR calibration is completed.
if (RUN_SALAMANDERS_NB2_ENV) all_dharma[["Salamanders_NB2_env"]] <- dharma_for_model(fit_env, "Salamanders", "NB2_env")
if (RUN_OWLS_ZINB_ARRIVAL) all_dharma[["Owls_ZINB_arrival"]] <- dharma_for_model(fit_arr, "Owls", "ZINB_arrival")
if (RUN_EPIL_NB2_REFINED) all_dharma[["epil_NB2_refined"]] <- dharma_for_model(fit_ref, "epil", "NB2_refined")

# ---------------- Final outputs ----------------
if (length(all_rr_summaries) > 0) write_csv(bind_rows(all_rr_summaries), file.path(OUT, "data",
    "RR_guided_refinement_summary_B1000.csv"))
if (length(all_dharma) > 0) write_csv(bind_rows(all_dharma), file.path(OUT, "data",
    "RR_guided_refinement_DHARMa_summary_B1000.csv"))
if (length(all_aic) > 0) write_csv(bind_rows(all_aic), file.path(OUT, "data",
    "RR_guided_refinement_AIC_summary.csv"))

sink(file.path(OUT, "logs", "sessionInfo_RR_guided_model_refinement_B1000.txt"))
print(sessionInfo())
sink()

cat("\nDONE.\n")
cat("All outputs are in:\n", OUT, "\n")
cat("\nOutput folder:\nRR_guided_model_refinement_B1000_submission_results\n")

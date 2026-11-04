# ============================================================
# Remaining empirical examples for GLMM recursive residual diagnostics
# B = 1000
#
# Examples included:
#   1. MASS::bacteria        binary longitudinal GLMM
#   2. glmmTMB::Owls         count / zero-inflated count GLMM
#   3. MASS::epil            epilepsy longitudinal count GLMM
#   4. datasets::ChickWeight Gamma positive-continuous GLMM
#   5. lme4::sleepstudy      Gamma positive-continuous GLMM
#   6. mlmRev::Contraception optional larger binary multilevel GLMM
#
# ============================================================

required <- c("lme4", "glmmTMB", "DHARMa", "statmod", "MASS",
              "dplyr", "ggplot2", "readr", "tibble")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install missing packages first: ", paste(missing, collapse = ", "))
}

library(lme4)
library(glmmTMB)
library(DHARMa)
library(statmod)
library(MASS)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

set.seed(2026)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

B_ENV <- 1000
N_QUAD <- 40
ALPHA <- 0.05

# Toggle examples. Set any to FALSE to skip.
RUN_BACTERIA      <- TRUE
RUN_OWLS          <- TRUE
RUN_EPIL          <- TRUE
RUN_CHICKWEIGHT   <- TRUE
RUN_SLEEPSTUDY    <- TRUE
RUN_CONTRACEPTION <- FALSE   # Optional; requires package mlmRev.

# Gamma RR note:
# The Gamma RR implementation assumes glmmTMB's sigma(fit) is the Gamma dispersion phi
# in Var(Y|b)=phi*mu^2, hence shape=1/phi and scale=mu*phi.
# This is the usual GLM Gamma mean-variance convention, but we will review results carefully.
RUN_GAMMA_RR <- TRUE

# Function definitions used in this script
source(file.path("code", "functions", "functions_05_empirical_additional_examples.R"), local = TRUE)

OUT <- file.path(get_output_root(), "additional_empirical_results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "logs"), showWarnings = FALSE, recursive = TRUE)

cat("\nOutputs will be written to:\n", OUT, "\n\n")

# ============================================================
# Utility functions
# ============================================================


# ============================================================
# Random-intercept Bernoulli recursive residuals for lme4::glmer
#
# current:
#   1. Uses an explicit fixed-effect formula for model.matrix().
#   2. Handles boundary random-effect variance safely.
#   3. Simulation loop is robust to failed refits/residualizations.
# ============================================================


# ============================================================
# Count RR for glmmTMB random-intercept models: Poisson/NB2/ZINB
# ============================================================


# ============================================================
# Gamma RR for glmmTMB random-intercept Gamma(log) models
# ============================================================


# ============================================================
# Run examples
# ============================================================

all_summaries <- list()
delayed_dharma <- list() # Comparator diagnostics run only after all RR simulations.
binary_dharma_queue <- new.env(parent = emptyenv())
binary_dharma_queue$tasks <- list()

# ------------------------------------------------------------
# 1. MASS::bacteria - binary longitudinal example
# ------------------------------------------------------------

if (RUN_BACTERIA) {
  cat("\nLoading MASS::bacteria...\n")
  data(bacteria, package = "MASS")
  bacteria2 <- bacteria

  # Robust binary conversion.
  bacteria2$ybin <- if (is.factor(bacteria2$y) || is.character(bacteria2$y)) {
    as.integer(as.character(bacteria2$y) %in% c("y", "yes", "Y", "Yes", "1", "TRUE", "present"))
  } else {
    as.integer(bacteria2$y)
  }

  bacteria2$ID <- factor(bacteria2$ID)
  bacteria2$week_num <- suppressWarnings(as.numeric(as.character(bacteria2$week)))
  if (any(is.na(bacteria2$week_num))) {
    bacteria2$week_num <- as.numeric(bacteria2$week)
  }

  # Keep a moderate and interpretable model.
  bacteria2$trt <- factor(bacteria2$trt)
  bacteria2$week_fac <- factor(bacteria2$week)

  all_summaries[["bacteria"]] <- run_binary_lme4_rr(
    dataset_name = "bacteria",
    data = bacteria2,
    formula = ybin ~ trt * week_fac + (1 | ID),
    fixed_formula = ~ trt * week_fac,
    id_var = "ID",
    order_var = "week_num",
    y_var = "ybin",
    index_vars = c("week_num", "pred_prob", "trt"),
    diagnostic_queue = binary_dharma_queue
  )
}

# ------------------------------------------------------------
# 2. glmmTMB::Owls - count / zero-inflated count example
# ------------------------------------------------------------

if (RUN_OWLS) {
  cat("\nLoading glmmTMB::Owls...\n")
  data(Owls, package = "glmmTMB")

  Owls2 <- Owls |>
    mutate(
      Nest = factor(Nest),
      FoodTreatment = factor(FoodTreatment),
      SexParent = factor(SexParent),
      ArrivalTime_num = suppressWarnings(as.numeric(ArrivalTime)),
      logBroodSize = log(BroodSize)
    )

  fit_owls_pois <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                             offset(logBroodSize) + (1 | Nest),
                           data = Owls2,
                           family = poisson)

  fit_owls_nb2 <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                            offset(logBroodSize) + (1 | Nest),
                          data = Owls2,
                          family = nbinom2)

  fit_owls_zip <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                            offset(logBroodSize) + (1 | Nest),
                          ziformula = ~ FoodTreatment * SexParent,
                          data = Owls2,
                          family = poisson)

  fit_owls_zinb <- glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                             offset(logBroodSize) + (1 | Nest),
                           ziformula = ~ FoodTreatment * SexParent,
                           data = Owls2,
                           family = nbinom2)

  write_model_AIC(
    list(Poisson_RI = fit_owls_pois, NB2_RI = fit_owls_nb2,
         ZIP_RI = fit_owls_zip, ZINB_RI = fit_owls_zinb),
    file.path(OUT, "data", "Owls_model_AIC.csv")
  )

  delayed_dharma[[length(delayed_dharma) + 1L]] <- function() {
  dharma_summary(
    list(Poisson_RI = fit_owls_pois, NB2_RI = fit_owls_nb2,
         ZIP_RI = fit_owls_zip, ZINB_RI = fit_owls_zinb),
    prefix = "Owls"
  )
  }

  owls_fixed <- ~ FoodTreatment * SexParent
  owls_zi <- ~ FoodTreatment * SexParent
  owls_offset <- function(d) d$logBroodSize
  owls_refit_pois <- function(d) glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                                           offset(logBroodSize) + (1 | Nest),
                                         data = d, family = poisson)
  owls_refit_nb2 <- function(d) glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                                          offset(logBroodSize) + (1 | Nest),
                                        data = d, family = nbinom2)
  owls_refit_zinb <- function(d) glmmTMB(SiblingNegotiation ~ FoodTreatment * SexParent +
                                           offset(logBroodSize) + (1 | Nest),
                                         ziformula = ~ FoodTreatment * SexParent,
                                         data = d, family = nbinom2)

  all_summaries[["Owls_Poisson"]] <- run_count_rr_model(
    "Owls", "Poisson_RI", fit_owls_pois, Owls2, owls_fixed,
    owls_refit_pois, "Nest", "ArrivalTime_num", "SiblingNegotiation",
    "poisson", FALSE, NULL, owls_offset,
    c("ArrivalTime_num", "pred_mean", "FoodTreatment", "SexParent", "BroodSize")
  )

  all_summaries[["Owls_NB2"]] <- run_count_rr_model(
    "Owls", "NB2_RI", fit_owls_nb2, Owls2, owls_fixed,
    owls_refit_nb2, "Nest", "ArrivalTime_num", "SiblingNegotiation",
    "nbinom2", FALSE, NULL, owls_offset,
    c("ArrivalTime_num", "pred_mean", "FoodTreatment", "SexParent", "BroodSize")
  )

  all_summaries[["Owls_ZINB"]] <- run_count_rr_model(
    "Owls", "ZINB_RI", fit_owls_zinb, Owls2, owls_fixed,
    owls_refit_zinb, "Nest", "ArrivalTime_num", "SiblingNegotiation",
    "nbinom2", TRUE, owls_zi, owls_offset,
    c("ArrivalTime_num", "pred_mean", "FoodTreatment", "SexParent", "BroodSize")
  )
}

# ------------------------------------------------------------
# 3. MASS::epil - epilepsy longitudinal count example
# ------------------------------------------------------------

if (RUN_EPIL) {
  cat("\nLoading MASS::epil...\n")
  data(epil, package = "MASS")
  epil2 <- epil

  if (!"subject" %in% names(epil2)) {
    if ("id" %in% names(epil2)) epil2$subject <- epil2$id
    if ("ID" %in% names(epil2)) epil2$subject <- epil2$ID
  }

  epil2$subject <- factor(epil2$subject)
  epil2$period_num <- suppressWarnings(as.numeric(as.character(epil2$period)))
  if (any(is.na(epil2$period_num))) epil2$period_num <- as.numeric(epil2$period)
  epil2$period_fac <- factor(epil2$period)
  epil2$trt <- factor(epil2$trt)

  fit_epil_pois <- glmmTMB(y ~ trt + base + age + period_fac + (1 | subject),
                           data = epil2,
                           family = poisson)

  fit_epil_nb2 <- glmmTMB(y ~ trt + base + age + period_fac + (1 | subject),
                          data = epil2,
                          family = nbinom2)

  write_model_AIC(
    list(Poisson_RI = fit_epil_pois, NB2_RI = fit_epil_nb2),
    file.path(OUT, "data", "epil_model_AIC.csv")
  )

  delayed_dharma[[length(delayed_dharma) + 1L]] <- function() {
  dharma_summary(
    list(Poisson_RI = fit_epil_pois, NB2_RI = fit_epil_nb2),
    prefix = "epil"
  )
  }

  epil_fixed <- ~ trt + base + age + period_fac
  epil_refit_pois <- function(d) glmmTMB(y ~ trt + base + age + period_fac + (1 | subject),
                                         data = d, family = poisson)
  epil_refit_nb2 <- function(d) glmmTMB(y ~ trt + base + age + period_fac + (1 | subject),
                                        data = d, family = nbinom2)

  all_summaries[["epil_Poisson"]] <- run_count_rr_model(
    "epil", "Poisson_RI", fit_epil_pois, epil2, epil_fixed,
    epil_refit_pois, "subject", "period_num", "y",
    "poisson", FALSE, NULL, NULL,
    c("period_num", "pred_mean", "base", "age", "trt")
  )

  all_summaries[["epil_NB2"]] <- run_count_rr_model(
    "epil", "NB2_RI", fit_epil_nb2, epil2, epil_fixed,
    epil_refit_nb2, "subject", "period_num", "y",
    "nbinom2", FALSE, NULL, NULL,
    c("period_num", "pred_mean", "base", "age", "trt")
  )
}

# ------------------------------------------------------------
# 4. datasets::ChickWeight - Gamma positive continuous example
# ------------------------------------------------------------

if (RUN_CHICKWEIGHT && RUN_GAMMA_RR) {
  cat("\nLoading datasets::ChickWeight...\n")
  data(ChickWeight, package = "datasets")
  chick <- ChickWeight |>
    mutate(Chick = factor(Chick),
           Diet = factor(Diet),
           Time_num = as.numeric(Time))

  fit_chick_gamma <- glmmTMB(weight ~ poly(Time_num, 2) * Diet + (1 | Chick),
                             data = chick,
                             family = Gamma(link = "log"))

  fit_chick_gauss <- lme4::lmer(log(weight) ~ poly(Time_num, 2) * Diet + (1 | Chick),
                                data = chick)

  aic_chick <- AIC(fit_chick_gamma, fit_chick_gauss) |>
    as.data.frame() |>
    tibble::rownames_to_column("model")
  readr::write_csv(aic_chick, file.path(OUT, "data", "ChickWeight_model_AIC.csv"))

  delayed_dharma[[length(delayed_dharma) + 1L]] <- function() {
  dharma_summary(list(Gamma_RI = fit_chick_gamma), prefix = "ChickWeight")
  }

  chick_fixed <- ~ poly(Time_num, 2) * Diet
  chick_refit <- function(d) glmmTMB(weight ~ poly(Time_num, 2) * Diet + (1 | Chick),
                                    data = d,
                                    family = Gamma(link = "log"))

  all_summaries[["ChickWeight_Gamma"]] <- run_gamma_rr_model(
    "ChickWeight", "Gamma_RI", fit_chick_gamma, chick, chick_fixed,
    chick_refit, "Chick", "Time_num", "weight",
    c("Time_num", "pred_mean", "Diet")
  )
}

# ------------------------------------------------------------
# 5. lme4::sleepstudy - Gamma positive continuous example
# ------------------------------------------------------------

if (RUN_SLEEPSTUDY && RUN_GAMMA_RR) {
  cat("\nLoading lme4::sleepstudy...\n")
  data(sleepstudy, package = "lme4")
  sleep <- sleepstudy |>
    mutate(Subject = factor(Subject),
           Days_num = as.numeric(Days))

  fit_sleep_gamma <- glmmTMB(Reaction ~ Days_num + (1 | Subject),
                             data = sleep,
                             family = Gamma(link = "log"))

  fit_sleep_gauss <- lme4::lmer(log(Reaction) ~ Days_num + (1 | Subject),
                                data = sleep)

  aic_sleep <- AIC(fit_sleep_gamma, fit_sleep_gauss) |>
    as.data.frame() |>
    tibble::rownames_to_column("model")
  readr::write_csv(aic_sleep, file.path(OUT, "data", "sleepstudy_model_AIC.csv"))

  delayed_dharma[[length(delayed_dharma) + 1L]] <- function() {
  dharma_summary(list(Gamma_RI = fit_sleep_gamma), prefix = "sleepstudy")
  }

  sleep_fixed <- ~ Days_num
  sleep_refit <- function(d) glmmTMB(Reaction ~ Days_num + (1 | Subject),
                                    data = d,
                                    family = Gamma(link = "log"))

  all_summaries[["sleepstudy_Gamma"]] <- run_gamma_rr_model(
    "sleepstudy", "Gamma_RI", fit_sleep_gamma, sleep, sleep_fixed,
    sleep_refit, "Subject", "Days_num", "Reaction",
    c("Days_num", "pred_mean")
  )
}

# ------------------------------------------------------------
# 6. Optional mlmRev::Contraception - larger binary multilevel example
# ------------------------------------------------------------

if (RUN_CONTRACEPTION) {
  if (!requireNamespace("mlmRev", quietly = TRUE)) {
    warning("RUN_CONTRACEPTION is TRUE, but package mlmRev is not installed. Skipping.")
  } else {
    cat("\nLoading mlmRev::Contraception...\n")
    data(Contraception, package = "mlmRev")
    contra <- mlmRev::Contraception

    # Inspect/adjust names if your installed data version differs.
    # Common variables include use, age, urban, livch, district.
    if (!"district" %in% names(contra)) stop("Cannot find district column in Contraception.")
    if (!"use" %in% names(contra)) stop("Cannot find use column in Contraception.")

    contra$district <- factor(contra$district)
    contra$row_order <- seq_len(nrow(contra))
    contra$usebin <- if (is.factor(contra$use) || is.character(contra$use)) {
      as.integer(as.character(contra$use) %in% c("Y", "Yes", "yes", "use", "1", "TRUE"))
    } else {
      as.integer(contra$use)
    }

    contra$urban <- factor(contra$urban)
    if ("livch" %in% names(contra)) contra$livch <- factor(contra$livch)

    formula_contra <- if ("livch" %in% names(contra)) {
      usebin ~ age + I(age^2) + urban + livch + (1 | district)
    } else {
      usebin ~ age + I(age^2) + urban + (1 | district)
    }

    fixed_formula_contra <- if ("livch" %in% names(contra)) {
      ~ age + I(age^2) + urban + livch
    } else {
      ~ age + I(age^2) + urban
    }

    all_summaries[["Contraception"]] <- run_binary_lme4_rr(
      dataset_name = "Contraception",
      data = contra,
      formula = formula_contra,
      fixed_formula = fixed_formula_contra,
      id_var = "district",
      order_var = "row_order",
      y_var = "usebin",
      index_vars = c("row_order", "pred_prob", "age", "urban"),
      diagnostic_queue = binary_dharma_queue
    )
  }
}

# Run DHARMa after all RR calibration simulations to prevent version-dependent
# comparator side effects from changing the RR reference distributions.
for (diagnostic in delayed_dharma) diagnostic()
for (diagnostic in binary_dharma_queue$tasks) diagnostic()

# ============================================================
# Final combined summary
# ============================================================

if (length(all_summaries) > 0) {
  final_summary <- dplyr::bind_rows(all_summaries)
  readr::write_csv(final_summary, file.path(OUT, "data",
                                            "ALL_remaining_examples_RR_summary_B1000.csv"))
}

sink(file.path(OUT, "logs", "sessionInfo_remaining_examples_B1000.txt"))
print(sessionInfo())
sink()

cat("\nDONE.\n")
cat("All outputs are in:\n")
cat(OUT, "\n")
cat("\nOutput folder:\n")
cat("additional_empirical_results\n")

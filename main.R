# Single entrypoint for the archived GLMM RR reproducibility package.
# Launch with Rscript main.R from the extracted ZIP root.
# Defaults to a safe validation-only mode. It never changes the working directory.
# Usage: Rscript main.R check
#        Rscript main.R smoke
#        Rscript main.R full
# IMPORTANT: full runs can require millions of model fits.
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) tolower(args[[1L]]) else "check"
if (!mode %in% c("check", "smoke", "full")) stop("Mode must be check, smoke or full")
required <- c("README.md", "code/01_simulation_bernoulli_power.R", "data/simulation_catalogue/simulation_catalogue_primary_B1000.csv")
missing <- required[!file.exists(required)]
if (length(missing)) stop("Run from extracted archive root. Missing: ", paste(missing, collapse=", "))
cat("GLMM RR archive preflight; mode:",mode,"\n")
if (mode == "check") {
  cat("Archive files present. No numerical experiments executed.\n")
  print(sessionInfo())
  quit(save="no", status=0)
}
if (mode == "smoke") {
  Sys.setenv(RR_GLMM_RUN_MODE="SMOKE_TEST",FAST_TEST="TRUE",RR_GLMM_OUTPUT_ROOT="results",OUT_ROOT="results")
  jobs <- c("code/01_simulation_bernoulli_power.R", "code/02_simulation_bernoulli_dharma_comparison.R")
  cat("SMOKE mode: tests only scripts 01 and 02; does not recreate paper's published numbers.\n")
} else {
  Sys.setenv(RR_GLMM_RUN_MODE="FULL_Q1",FAST_TEST="FALSE",RR_GLMM_OUTPUT_ROOT="results",OUT_ROOT="results")
  jobs <- sprintf("code/%02d_%s.R",c(1,2,3,4,5,6,8,9),c("simulation_bernoulli_power","simulation_bernoulli_dharma_comparison","empirical_epilepsy_refinement","empirical_salamanders_analysis","empirical_additional_examples","ablation_epilepsy","simulation_catalogue","empirical_cbpp"))
}
for (job in jobs) {
  cat("\n== Running ",job," ==\n",sep="")
  t0 <- Sys.time()
  tryCatch(source(job,local=new.env(parent=globalenv())),error=function(e) stop("Failure in ",job,": ",conditionMessage(e)))
  cat("Completed in",as.numeric(difftime(Sys.time(),t0,units="mins")),"minutes\n")
  if (mode == "smoke") {
    if (grepl("01_simulation", job, fixed=TRUE)) {
      path <- "results/bernoulli_power_results/data/bernoulli_Q1_replicate_results.csv"
      dat <- read.csv(path)
      if (!nrow(dat) || any(!is.finite(dat$p_combo_minP)) || any(is.na(dat$detected_combo)))
        stop("SMOKE NUMERICAL VALIDATION FAILED: Bernoulli RR calibration has missing/invalid p-values; inspect ",path)
      cat("RR numeric gate: Bernoulli ",nrow(dat),"/",nrow(dat)," finite p-values\n",sep="")
    }
    if (grepl("02_simulation", job, fixed=TRUE)) {
      path <- "results/bernoulli_dharma_comparison_results/data/dharma_head_to_head_replicate_results.csv"
      dat <- read.csv(path)
      if (!nrow(dat) || any(!is.finite(dat$rr_pcombo)) || any(is.na(dat$rr_detect)))
        stop("SMOKE NUMERICAL VALIDATION FAILED: DHARMa comparison has missing/invalid RR statistics; inspect ",path)
      tests <- cbind(dat$dharma_uniformity_p, dat$dharma_dispersion_p, dat$dharma_outlier_p)
      if (any(!is.finite(tests)) || any(tests < 0 | tests > 1) ||
          any(abs(dat$dharma_min_raw_p - apply(tests, 1, min)) > 1e-12) ||
          any(abs(dat$dharma_min_bonf_p - pmin(1, 3 * apply(tests, 1, min))) > 1e-12))
        stop("SMOKE DHARMa THREE-TEST VALIDATION FAILED: inspect ",path)
      cat("DHARMa numeric gate: ",nrow(dat),"/",nrow(dat)," valid three-test decisions\n",sep="")
      cat("RR numeric gate: comparison ",nrow(dat),"/",nrow(dat)," finite RR p-values\n",sep="")
    }
  }
}
if (mode == "full") cat("Optional table formatter: python3 code/07_format_tables.py\n")

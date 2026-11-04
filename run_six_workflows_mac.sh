#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(dirname "$0")"
MODE="${1:-quick}"
mkdir -p workflow_verification/logs workflow_verification/runs
case "$MODE" in
 quick) JOBS=(08 09) ;;
 empirical) JOBS=(03 04 05 06) ;;
 all) JOBS=(08 09 03 04 05 06) ;;
 03|04|05|06|08|09) JOBS=("$MODE") ;;
 *) echo "Usage: ./run_six_workflows_mac.sh [quick|empirical|all|03|04|05|06|08|09]"; exit 2 ;;
esac
# Fail fast on syntax defects before copying data or starting any computations.
Rscript -e 'f <- sort(list.files("code", full.names=TRUE)); f <- f[endsWith(f,".R") & grepl("^[0-9][0-9]_",basename(f))]; stopifnot(length(f)>=8L); for (x in f) { tryCatch(parse(file=x), error=function(e) stop("R_PARSE_FAILED ",x,": ",conditionMessage(e))); cat("PARSE_OK ",x,"\n",sep="") }' | tee workflow_verification/logs/parse_preflight.log
Rscript -e 'pkgs<-c("lme4","glmmTMB","DHARMa","statmod","MASS","dplyr","ggplot2","readr","tibble"); m<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop("Missing: ",paste(m,collapse=",")); cat("DEPENDENCIES_OK\\n")' | tee workflow_verification/logs/dependencies.log
Rscript -e 'pkgs<-c("lme4","glmmTMB","DHARMa","statmod","MASS","dplyr","ggplot2","readr","tibble");invisible(lapply(pkgs,require,character.only=TRUE));writeLines(capture.output(sessionInfo()),"workflow_verification/logs/sessionInfo_loaded.txt")'
for ID in "${JOBS[@]}"; do
 case "$ID" in
 03) F=03_empirical_epilepsy_refinement.R ;;
 04) F=04_empirical_salamanders_analysis.R ;;
 05) F=05_empirical_additional_examples.R ;;
 06) F=06_ablation_epilepsy.R ;;
 08) F=08_simulation_catalogue.R ;;
 09) F=09_empirical_cbpp.R ;;
 esac
 WORK="workflow_verification/runs/$ID"
 mkdir -p "$WORK/results"
 for D in code data tables; do
  if [[ -L "$WORK/$D" ]]; then rm "$WORK/$D"; fi
  if [[ ! -d "$WORK/$D" ]]; then cp -R "$D" "$WORK/$D"; fi
 done
 echo "START $ID $(date -u +%FT%TZ)" | tee -a workflow_verification/logs/status.tsv
 START=$SECONDS
 (cd "$WORK" && Rscript "code/$F") > "workflow_verification/logs/$ID.log" 2>&1 && RC=0 || RC=$?
 echo -e "$ID\t$RC\t$((SECONDS-START))\t$(find "$WORK/results" -type f | wc -l | tr -d ' ')" | tee -a workflow_verification/logs/summary.tsv
 if [[ "$RC" -ne 0 ]]; then echo "FAILED $ID: inspect workflow_verification/logs/$ID.log"; exit "$RC"; fi
 if [[ "$(find "$WORK/results" -type f | wc -l)" -eq 0 ]]; then echo "FAILED $ID: zero result files"; exit 3; fi
 echo "PASS $ID ($((SECONDS-START)) seconds)"
done
python3 verify_six_workflows.py

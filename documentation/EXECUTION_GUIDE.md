# Execution instructions

**Starting location:** extract the single ZIP and run commands from `GLMM_Recursive_Residuals_Code_and_Data/`. Do not use absolute host-specific paths. The article itself is not included; all source files are provided for computational use.

## Verify supplied computational results without Monte Carlo resimulation

```bash
cd GLMM_Recursive_Residuals_Code_and_Data
python3 verify_results.py
Rscript main.R check
python3 code/07_format_tables.py  # WARNING: overwrites output .tex tables; do this in a working COPY only
```

`verify_results.py` performs the archived-versus-latest-targeted 76-file CSV comparison; this does not execute the original R experiments. `main.R check` checks only selected file existence and displays `sessionInfo()`, not scientific numerical reproducibility. `code/07_format_tables.py` is a formatter, not a statistical reanalysis. Run only inside an extracted working **copy** if preserving distributed file hashes is necessary.

## Limited functionality check

```bash
Rscript main.R smoke
```

Runs scripts 01 and 02 with reduced simulation settings, not manuscript-strength output. It may take noticeable time. Do **not** equate a smoke-test pass with exact reproduction of published results. The generated `results/` tree is excluded from the archive; work in a copy if desired.

## Complete analysis order (computationally expensive)

1. `Rscript code/01_simulation_bernoulli_power.R` — Bernoulli power/null-level simulation.
2. `Rscript code/02_simulation_bernoulli_dharma_comparison.R` — head-to-head RR/DHARMa simulation.
3. `Rscript code/03_empirical_epilepsy_refinement.R` — Salamanders/Owls/epilepsy refinement, B=1000.
4. `Rscript code/04_empirical_salamanders_analysis.R` — Salamanders Poisson/NB2/ZINB, B=1000.
5. `Rscript code/05_empirical_additional_examples.R` — additional empirical datasets, B=1000.
6. `Rscript code/06_ablation_epilepsy.R` — epilepsy ablation, Part A + costly Part B.
7. `python3 code/07_format_tables.py` — generate manuscript-oriented LaTeX ablation/calibration tables from recorded CSVs.
8. `Rscript code/08_simulation_catalogue.R` — format/verify **archived** simulation catalogue outputs; NOT a rerun of their scenario-level generators.
9. `Rscript code/09_empirical_cbpp.R` — format/verify **archived** CBPP RR/DHARMa results; NOT a fresh CBPP B=1000 calibration.

An alternate controller is `Rscript main.R full` (scripts 01,02,03,04,05,06,08,09 sequentially; run the Python formatter separately). `run_six_workflows_mac.sh` is an optional controller for scripts 03,04,05,06,08,09 and is **not** a complete equivalent of the nine-step sequence. Both full modes are computationally expensive; use a separate writable working copy. Note script 06 has `SMOKE_TEST <- FALSE` internally and requires manually setting that flag to run a short ablation check.

## Dependencies and built-in datasets

R scripts use base R and CRAN packages `lme4`, `glmmTMB`, `DHARMa`, `statmod`, `MASS`, `dplyr`, `ggplot2`, `readr`, `tibble`; optional `mlmRev` is used only if available in the additional-examples script. Python 3 standard library is used for table formatting and verification. The original simulations were run on a Mac with the historical `DHARMa 0.4.7` for empirical diagnostic comparisons. The precise previously recorded versions include `glmmTMB 1.1.14`, `lme4 2.0.1`, `DHARMa 0.4.7`, `MASS 7.3.65`, `statmod 1.5.2`, `ggplot2 4.0.3`, `dplyr 1.2.1`, `readr 2.2.0`, `tibble 3.3.1`. See `verification/logs/historical_sessionInfo.txt` for recorded environment and run `Rscript -e 'sessionInfo()'` on the reader's own system. Exact numerical compatibility can depend on R, package and numerical-library versions.

Observed datasets are taken from installed R packages: `glmmTMB::Salamanders`, `glmmTMB::Owls`, `MASS::epil`, `MASS::bacteria`, `datasets::ChickWeight`, `lme4::sleepstudy`, and optionally `mlmRev::Contraception`; CBPP is available from `lme4`. `data/` contains archived computed inputs/reference outputs, including the CBPP and simulation catalogue tables. The code generating **all** those archived source simulations has not been independently established in this particular supplement; scripts 08 and 09 reproduce derived tables from saved inputs.

## Numerical comparison evidence and scope

- In the full-scale empirical comparison, **69/69 RR/other files** matched archived numerical/text values; zero differing cells among **267,179 numerical cells** at abs tolerance `1e-7` (68 byte-identical).
- Targeted reproductions of DHARMa summaries for workflows 03/04 and 3/5 workflow-05 diagnostic files gave **three** remaining bootstrap p-value discrepancies: Owls NB2 `1.00 → 0.86`, epilepsy Poisson `0.22 → 0.24`, epilepsy NB2 `0.52 → 0.56`, all unchanged at alpha=0.05. Original reference p-values are retained.
- The corrected workflow-03 reporting path passed targeted execution after an R offset-function name collision in the original full-scale controller. A fresh **uninterrupted full** B=1000 replay of corrected script 03 was not performed.
- The finalized main paper's **14 included table files have been byte-compared** to the original final manuscript source and match. Its three numbered PDF figures are included under `paper_figures/`; a complete manuscript-to-every-numeric-claim audit has not been performed.

## Local file paths

Run from the extracted root; R scripts generally use root-relative paths or construct paths from the current working directory. Do not execute in a different working directory. The optional `run_six_workflows_mac.sh` locates its own directory using `cd` and then uses relative paths internally; this is documented rather than concealed.

# Code and data accompanying “Posterior-Predictive Recursive Residuals for GLMMs: Construction, Theory, and Calibrated Diagnostics”

**Research article:** *Biometrical Journal* · manuscript identifier 1725547  
**Contents:** R and Python analysis code, datasets and saved computational results, manuscript-linked tables and figures, a prototype R package, and numerical reproducibility evidence.

## Scope

This archive allows researchers to inspect the data and computational methods, reproduce the supplied analyses where source-generating code is available, and independently examine the numerical results. The published manuscript is not embedded in this archive. Neither the manuscript nor its archived scientific outputs have been revised for this release.

The materials distinguish **independent computation from source-level regeneration**. In particular, the saved 41-row simulation catalogue can be examined and its reporting tables regenerated, but its complete original scenario-generating program and data-generating specifications are not included. A text description of 40 scenarios in the article and 41 rows in the saved catalogue should not be assumed to be equivalent without identifying the additional scenario. The CBPP calibration was independently reconstructed; its overall calibrated inference agrees with the archived result, while a small number of statistics from singular fitted models are not identical to the historical samples.

## Directory guide

| Location | Contents |
|---|---|
| `code/` | Numbered R programs, helper functions, and a Python table formatter |
| `data/` | Saved empirical results, calibration samples and simulation-summary inputs |
| `tables/` | LaTeX tables related to the article and supplementary analyses |
| `figures/`, `paper_figures/` | Graphical outputs and article-linked figures |
| `R_package/` | `glmmRRdiag` source, documentation, unit tests and package archive |
| `documentation/` | Execution order, dependencies, hardware/software environment, random seeds, runtime and table/figure crosswalk |
| `verification/` | Original and regenerated results, comparison files, CBPP boundary analysis and detailed computational logs |
| `main.R`, `run_six_analysis programs_mac.sh` | Analysis entry points |
| `verify_results.py`, `verify_archived_results.py`, `verify_numerical_reproducibility.py` | Numerical validation utilities |
| `SHA256SUMS.txt` | SHA-256 integrity hashes for every other file |

See `FILE_GUIDE.md` and `documentation/EXECUTION_GUIDE.md` for the location and execution purpose of individual files.

## Getting started

R and Python 3 are required for computational checks. The R analyses use the CRAN packages `lme4`, `glmmTMB`, `DHARMa`, `statmod`, `MASS`, `dplyr`, `ggplot2`, `readr` and `tibble`; some datasets are distributed by those R packages rather than duplicated in `data/`. Further version details appear in `documentation/EXECUTION_GUIDE.md`, `documentation/R_SESSION_INFO.txt`, and the recorded sessions under `verification/logs/`.

Extract the archive into any writable folder and execute the commands **from its root directory**:

```bash
cd GLMM_Recursive_Residuals_Code_and_Data
shasum -a 256 -c SHA256SUMS.txt
python3 verify_results.py
Rscript main.R check
```

On systems lacking macOS `shasum`, Python can verify the same hashes using `hashlib`. The checksum command is a file-integrity test, while `verify_results.py` compares saved numerical outputs. `Rscript main.R check` checks selected inputs and prints the R session; it does not rerun all simulations.

For a limited functionality check:

```bash
Rscript main.R smoke
```

For the full analysis sequence and estimated computing burden, use `documentation/EXECUTION_GUIDE.md`. Run computationally expensive scripts **in a copy**: they can create or overwrite files, and rerunning all numerical analyses is not necessary merely to inspect the deposited results. `code/08_simulation_catalogue.R` and `code/09_empirical_cbpp.R` format saved results; these scripts are not the upstream generators for the corresponding calibration experiments.

## Reproduction evidence and interpretation

In the available empirical model-calibration evidence, 69 paired RR/related output files agree numerically, including 267,179 compared numerical entries at absolute tolerance 1e-7. Three DHARMa bootstrap outlier p-values differ (Owls NB2: 1.00 versus 0.86; epilepsy Poisson: 0.22 versus 0.24; epilepsy NB2: 0.52 versus 0.56), without changing the corresponding decisions at significance level 0.05.

The independent CBPP reconstruction produced 1,000 calibration refits with no fitting failures, 16 singular refits and a calibrated combined p-value of approximately 0.08591409, agreeing with the archived combined result. Of 10,000 compared replicate-statistic entries, 28 entries across 15 singular refits were historically nonidentical; the other eight statistic columns matched numerically. A mathematically justified boundary convention was evaluated for zero fitted random-effect variance. It preserves the overall calibrated result but **is not a recovery of the exact historical numerical ordering**.

A corrected name collision in the epilepsy reporting script is documented in the corresponding computational evidence. The corrected analysis path was checked using saved simulation outputs; a separate uninterrupted full B=1000 execution of that script was not performed. The origin of the residual DHARMa bootstrap differences has not been established. These qualifications are essential for interpreting numerical reproduction claims.

## Data use and scientific interpretation

The accompanying R code, result files, and documentation should be interpreted together with the published methods, model specifications, statistical definitions and figure/table captions. Manuscript-linked output paths are provided in `documentation/TABLE_FIGURE_CROSSWALK.csv`. The archived results are preserved as reference values, not silently replaced with outputs from reconstructed methods. For independent analyses, researchers should report software versions, parameter choices, random-number initialization, and any deviations from the archived models.

This archive does not imply that every experiment can be regenerated from raw simulated observations using only its contents. It provides executable analyses where available and transparent saved inputs and numerical comparisons elsewhere.

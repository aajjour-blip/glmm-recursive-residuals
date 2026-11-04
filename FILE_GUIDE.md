# File guide

The layout follows standard computational-research conventions. Paths below are relative to the extracted archive root.

| Location | Research purpose |
|---|---|
| `README.md` | Methods context, software requirements, instructions, interpretation |
| `START_HERE.md` | Short orientation for readers |
| `REPRODUCIBILITY.md` | Numerical agreement, limitations and model-specific evidence |
| `documentation/EXECUTION_GUIDE.md` | Script-by-script execution order, package dependencies, data sources |
| `documentation/COMPUTATIONAL_HARDWARE.txt` | Recorded computing platform |
| `documentation/R_SESSION_INFO.txt` | R environment and numerical libraries |
| `documentation/RUNTIME_AND_PARALLELIZATION.md` | Computing requirements and parallel execution |
| `documentation/RANDOM_SEEDS_AND_INTERMEDIATES.md` | RNG and saved intermediate-data information |
| `documentation/TABLE_FIGURE_CROSSWALK.csv` | Association of tables and figures with supporting outputs |
| `code/` | Research program source and helper functions |
| `data/` | Empirical reference data, stored simulation/calibration outputs |
| `tables/`, `figures/`, `paper_figures/` | Research tables, figures and exported plots |
| `R_package/glmmRRdiag/` | Prototype analysis package, tests and usage documentation |
| `verification/` | Paired results, computational logs, CBPP numerical studies and archived comparisons |
| `main.R` | R entry point for checking or running source programs |
| `run_six_workflows_mac.sh` | Optional Mac-compatible R controller for supported program groups |
| `verify_results.py` | Numerical agreement checks using stored reference and regenerated results |
| `SHA256SUMS.txt` | Hashes for archive integrity verification |

Some files inside `verification/` retain the filenames assigned by their original computation. They are preserved to maintain exact correspondence with generated outputs, command logs and numerical comparison data. They are supporting evidence, not separate versions of the article.

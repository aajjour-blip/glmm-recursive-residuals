# Research materials — reader's guide

These files accompany the research article **“Posterior-Predictive Recursive Residuals for GLMMs: Construction, Theory, and Calibrated Diagnostics.”** They are intended for researchers who wish to understand the computational methods, inspect underlying archived data and outputs, or reproduce the available numerical analyses.

**Start with `README.md`** for scope, dependencies, computational limitations and an interpretation of the saved evidence. **Use `FILE_GUIDE.md`** to locate code, datasets, tables, figures and model-calibration records. **Use `documentation/EXECUTION_GUIDE.md`** for commands and the required software environment.

For a quick read-only check after extracting the archive, run:

```bash
cd GLMM_Recursive_Residuals_Code_and_Data
shasum -a 256 -c SHA256SUMS.txt
python3 verify_results.py
```

Run the analysis scripts in a working copy if new simulations or regenerated tables are desired. The published scientific outputs remain available in the `data/`, `tables/` and `figures/` folders.

**Scope of reproduction:** The empirical RR comparisons and CBPP calibration conclusions are numerically supported by the included records. Some historical bootstrap statistics differ, and the upstream generator for the full simulation catalogue is not included. See `REPRODUCIBILITY.md` for the precise boundaries of independent reproduction.

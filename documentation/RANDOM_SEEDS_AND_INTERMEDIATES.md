# Seeds, checkpoint files, and fast reproducibility checks

The R sources contain fixed master seeds (usually `2026`) and per-replicate/per-calibration seeds. The Bernoulli simulation functions explicitly use identifiers to derive separate seeds. Workflow 06 Part B uses `base_seed = 1e7 * scenario_id + 1e4 * replicate_id` with separate offsets for simulated datasets, observations, calibration and reference distributions, enabling scenario/replicate-level deterministic replay given the same software implementation.

Workflows 03–05 produce **14 B=1000 simulation-statistic CSVs** with 1,000 replicates each. Their archived calibration arrays, component p-values and model summaries are included under `data/refinement/`, `data/empirical/`, `data/remaining/` and `verification/regenerated_output/`. The 69-file numerical agreement report allows independent spot checks without re-running full R jobs. Workflow 06 Part B intermediate scenario results are stored under `data/ablation_epilepsy/partB_*.rds`; Part A has `partA_results.rds`. Bernoulli and DHARMa head-to-head simulation replicate-level results are stored under `data/completed_simulations/`.

**Important qualification:** A fixed seed does not by itself guarantee identical bootstrap diagnostics across R/DHARMa versions or distinct RNG states. Three workflow-05 DHARMa bootstrap p-values are not exactly reproduced. Results should not be overwritten to make them pass. The 14 model calibrations and source-dependent empirical RR output files are numerically reproducible as documented.

Quick inspection (no long simulations):

```bash
python3 verify_results.py
Rscript main.R check
```

Smoke only for Bernoulli experiments: `Rscript main.R smoke` (not a full replacement). The original `code/01_...` and `code/02_...` scripts provide a small-run mode; for script 06 use its explicitly documented `SMOKE_TEST` flag in a working copy.

**Limitations:** `code/08_simulation_catalogue.R` consumes archived scenario-level results, and `code/09_empirical_cbpp.R` consumes archived CBPP summary/calibration CSVs. The source generators for these two sets of archived inputs are not demonstrated here; therefore do not claim ground-up reproduction of every underlying computation from scratch. Their derived manuscript tables are reproducible from saved intermediates.

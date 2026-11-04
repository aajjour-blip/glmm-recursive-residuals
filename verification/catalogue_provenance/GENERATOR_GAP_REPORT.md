# Simulation catalogue generator status

All 41 archived scenarios pass internal validation (16 Bernoulli, 21 Poisson, 4 Gamma). The original upstream generator is not established. Workflow 08 regenerates tables, not raw scenario calibrations.

Each scenario requires its original DGP, fitted working GLMM, cluster dimensions, covariate distributions, parameter values, simulation ordering, RNG seed schedule, recursive residual definition, calibration/refit failure treatment and omnibus calibration. These are **not recoverable uniquely from archived p-values**.

First resolve the NOT_ESTABLISHED fields in `SCENARIO_SPECIFICATION_MATRIX.csv` from primary code or detailed methods; then implement independently labelled reconstructed generators and validate reduced-run null/alternative pilots before full B=1000. Do not change reference results. The journal checklist all-code requirement remains open for ground-up catalogue generation.

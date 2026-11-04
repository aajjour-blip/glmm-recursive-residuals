# CBPP full-scale independent numerical comparison

Manuscript 1725547 — 4 November 2026

## Results

- Calibration evidence archive: integrity verified (13 files).
- Reference: archived CBPP empirical CSV in the journal package.
- 1000 rows × 10 statistic columns = 10,000 numerical cells compared with absolute tolerance 1e-7.
- 9,972/10,000 matching, **28 discrepancies affecting 15 replicates**.
- All 15 affected replicates are among the 16 reconstructed singular fits; singular replicate 831 has zero discrepancies.
- Affected columns: `Tinf_pred_prob` (13 cells; max absolute difference 0.598647163), `T2_pred_prob` (15 cells; max absolute difference 0.19594512).
- All eight other statistic columns have zero differences exceeding 1e-7.
- All 10 component calibration rows agree numerically, maximum difference 2e-15.
- Combined p-value: archived 0.08591408591408592; reconstructed 0.0859140859140859, non-rejection at 5% unchanged.
- Observed RR residuals: row-aligned, numerical matches within 1e-7 (inspected).
- 16 singular refits, zero failures (Mac checkpoint).

## Interpretation

The independently reconstructed CBPP generator reproduces the observed statistics, component calibration summaries and overall conclusion. It is **not replicate-wise identical** on singular fits. The exact reason the predictive-probability components diverge requires inspection of the singular-fit numerical branch and the generator's original historical implementation. This is consistent with, but does not prove, differences in singularity handling. Do not claim 10,000/10,000 exact numerical agreement.

The original upstream CBPP generator was not recovered. The reconstructed implementation is a new, validated candidate and should be identified as such. The 41-scenario simulation catalogue is separately unresolved.

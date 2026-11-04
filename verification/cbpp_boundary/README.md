# CBPP singular-boundary resolution (manuscript 1725547)

This is a **prospectively defined numerical correction**, not a recovery of the exact historical generator. Preserve the finalized manuscript and all archived reference files unchanged.

## Why

In 15 of the 16 singular CBPP refits, estimated random-intercept SD is exactly zero. Quadrature-based posterior predictive probability calculations produce 5–7 distinct floating-point values despite the mathematical prediction being identical within a period; these near-ties affect `Tinf_pred_prob` and `T2_pred_prob`.

## Rule

Only for `sigma == 0` (exact boundary), compute the *index* probability using the analytically degenerate Gaussian random-intercept integral, `plogis(X %*% fixef(fit))`. Use stable ascending sorting of this index, then original observation order for exact ties. The randomized PIT residuals and all other components remain UNMODIFIED. For near-zero but nonzero SD, preserve original algorithm (including singular replicate 831). No historical values are overwritten. Apply this same rule to every refit (and observed model if its sigma is exactly zero) in any future full prospective run.

The combined analysis in this package corrects only the 15 saved singular refits, leaving 985 other simulations unchanged, then recalibrates the entire 1,000-replicate mixture using the package's original `calibrate_components()` function. This is valid for the explicit boundary-branch algorithm but it must be described as a **new convention**, not historical replicate reproduction.

## Scope and caveats

* Full-precision R reconstruction of the 16 singular replicates verifies the original previous checkpoint first (tolerance 1e-7) and the other eight statistics, then changes only the index of the two probability-based statistics in exactly 15 zero-variance refits.
* A separate Python CSV-precision preview cross-checks the direction of corrections; the exported CSV truncates floating-point information and cannot certify historical ordering.
* It does **not** establish the ordering used by the unavailable historical 41-scenario generator, or make archived 28 discrepancies disappear.
* The original archived CBPP p-value and manuscript stay unchanged. The new p-value is separately labeled and statistically interpreted only after your Mac run.

## Mac command

```bash
cd ~/Downloads && \
unzip -qo GLMM_1725547_CBPP_BOUNDARY_RESOLUTION.zip && \
cd GLMM_1725547_CBPP_BOUNDARY_RESOLUTION && \
chmod +x run_mac.sh && \
./run_mac.sh
```

Output: `GLMM_1725547_CBPP_BOUNDARY_RESOLUTION_RESULTS.zip`. The R script requires `lme4`, `statmod`, and `dplyr`.

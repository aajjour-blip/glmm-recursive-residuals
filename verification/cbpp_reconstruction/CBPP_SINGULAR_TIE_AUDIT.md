# Independent CBPP singular-refit numerical audit — 4 November 2026

## Inputs and scope

The six saved diagnostic records were examined, including row-level recursive residuals, predicted probabilities, fit diagnostics, and tested ordering conventions. Cross-checked reconstructed-versus-archived discrepancies against the source program. No models were refitted in this independent audit.

## Confirmed

- The 16 singular refits are individually reproducible by the Mac targeted runner to numerical precision approximately 5e-15 for their recorded reconstructed statistics. Fifteen have estimated random-intercept standard deviation zero; replicate 831 is approximately 2e-7.
- Of 10,000 CBPP calibration statistics, 28 differences across 15 singular refits remain in `Tinf_pred_prob` (13) and `T2_pred_prob` (15). The other eight statistics match within 1e-7. Archived calibration p-values and omnibus decision match.
- With sigma zero, fitted predictive probabilities should be constant within each CBPP period under the fitted fixed-effect specification; the **in-memory diagnostic** reports 5–7 distinct exact predicted probabilities for each of the 15 zero-variance fits, versus exactly four values after rounding to 12 decimal places (one per period). This is direct evidence of floating-point deviations from exact ties, presumably generated in quadrature normalization or summation; the specific source operation has not yet been isolated. The saved row-level CSV serializes to approximately 15 significant digits and typically erases these very small within-period differences, so the raw CSV alone is **not** a lossless record of the actual in-memory ordering. `order(pred_prob, row)` can use those tiny differences to select different within-period permutations.
- The existing deterministic rounding and bucketing rules did not reproduce both archived probability-indexed statistics in any of the 15 mismatching refits. The one match was replicate 831, which matched before adjustments. It is therefore incorrect to say an alternative deterministic tie rule resolved a mismatch.

## Independent sensitivity diagnostic

For each of 16 replicates, drew 10,000 fixed-seed random permutations within the four CBPP period groups, preserving each residual and period (no refitting). Recalculated max absolute cumulative RR and mean squared cumulative RR. The output CSV records ranges, nearest joint distance to the archived pair, and counts with joint maximum absolute error under 0.01. Each archived statistic fell within the separately simulated range of its statistic, and at least one permutation had joint error below 0.01 for every replicate. This makes intra-period ordering **plausible** as the explanation; it does not reproduce the exact historical ordering or show exact joint attainability. **Important:** these random permutations are performed on rounded CSV exports of the residuals, not on the original full-precision prediction vectors.

Because the ordering of almost tied probabilities is implementation-sensitive, even a minuscule numerical perturbation can give a materially different cumulative-process statistic. The evidence implicates numerical ordering/tie sensitivity at the boundary rather than a change in the fitted mean or the randomized residual values. It does **not** establish that conditional/marginal prediction, quadrature normalization or the historical R library version can be excluded as contributing factors.

## Decision and remedy

Do not overwrite historical replicate statistics, tune a different tie-breaking rule to match them, or call the original generator exactly recovered. In a future *prospectively defined* implementation, use a documented, invariant tie policy (e.g. group probabilities under a declared numerical tolerance and break ties in fixed observation order), then regenerate both calibration and observed statistics consistently. That would be a **new convention**, not evidence of reproducing the archived historical stream. The finalized manuscript and archival numbers stay unchanged. The simulation-catalogue generator remains a separate unresolved item.

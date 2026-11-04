# Reproducibility and numerical interpretation

This document states what can be independently checked using the supplied research materials. It is intended for statistical researchers and users of the code and data.

## Empirical models

The archive contains the source programs, recorded inputs, saved B=1000 calibration results, numerical comparison files and corresponding figures and tables. Across the audited empirical RR results, 69/69 files agree within the documented numerical tolerance, including 267,179 compared numerical entries. The three remaining DHARMa bootstrap p-value differences are documented in `verification/dharma_bootstrap_differences.csv`; the significance classifications at the 5% threshold are unchanged. Some calibration programs are expensive and depend on software versions.

## CBPP grouped-binomial calibration

The independently reconstructed calibration includes 1,000 refits, 16 singular fits and no fit failures. The calibrated combined p-value agrees with the archived result (approximately 0.08591409). Of the 10,000 compared scalar replicate-statistic values, 28 values from 15 singular refits differ from the historical archive; two statistics indexed by predictive probabilities are affected. A zero-variance probability convention was evaluated and retains the combined calibration result. This is numerical sensitivity analysis, not evidence of exact reconstruction of the original historical simulation code. See `verification/cbpp_reconstruction/` and `verification/cbpp_boundary/`.

## Simulation catalogue

The saved simulation catalogue has 41 entries: 16 Bernoulli, 21 Poisson, and four Gamma scenarios. Code `code/08_simulation_catalogue.R` derives reporting tables from the saved catalogue inputs. It does **not** create the primary scenario-level simulated datasets. The complete historical scenario generator and parameter specifications were not recovered; its primary results therefore should be treated as archived inputs rather than independently regenerated simulation outcomes. The source article describes 40 scenarios; the relationship between that description and the 41-row archive remains unresolved. See `verification/catalogue_provenance/`.

## Practical use

For numerical verification without resimulation, execute `python3 verify_results.py` from the extracted archive root. To reproduce analysis outputs, consult `documentation/EXECUTION_GUIDE.md` and the individual script comments. Exact agreement of stochastic comparisons can depend on the software and RNG environment. Preserve the archived reference outputs when performing new calculations.

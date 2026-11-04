# Empirical numerical comparison


## Verification scope
- 76 CSV pairs from workflows 03, 04, and 05, archived versus regenerated.
- Absolute numeric tolerance: 1e-07.
- RR/other: **69/69 match**, 68 byte identical, 267,179 numerical cells, 0 numerical differences, 0 text differences.
- DHARMa: **2/7 match**, 5 differing files, 11 numerical discrepancies, 0 significance-decision changes at alpha=0.05.

## Caveats
- The original workflow-03 controller run failed during Owls DHARMa because of an `offset()` helper name collision; its RR results and reporting outputs were recovered. Do not describe it as an uninterrupted original run.
- 69/69 RR/other matching does not mean that all 76 files are exact, nor that all manuscript figures or theory have been independently validated.
- The 11 DHARMa diagnostic differences remain genuine and must not be silently changed in the archive.
- The archived reference outputs are retained unchanged.

## Contents
- `archived_reference/`: all 76 archived comparison CSVs.
- `regenerated_output/`: all 76 Mac-regenerated comparison CSVs.
- `complete_workflow_outputs/`: full Mac-generated output trees for workflows 03, 04, and 05, including ancillary files and figures.
- `audit/`: fresh verifier output and full 11-row DHARMa discrepancy report.
- `logs/`: available original controller logs (not evidence of a successful workflow-03 controller completion).
- `SHA256_MANIFEST.csv`: independent SHA-256 hashes for every archived and regenerated CSV included in the package.

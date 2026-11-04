# Runtime, parallelization, and hardware

## Empirical observed wall times

The recorded full-scale Mac controller ran the three independent empirical workflows concurrently (one process per workflow: **03, 04, 05**). It is **not** six-way parallel simulation. Individual successful workflow elapsed times reported by the Mac controller:

| Workflow | Recorded wall time | Detail |
|---|---:|---|
| 03 | Not available as an uninterrupted successful run | Full RR B=1000 calculations finished; original run stopped during DHARMa; completed by targeted reporting recovery |
| 04 | 1,583 s (about 26 min) | Three empirical RR B=1000 models |
| 05 | 2,853 s (about 48 min) | Eight empirical RR B=1000 models |
| 06 | Approximately 35,316 s (about 9.8 h), from an earlier Mac execution | Ablation simulation: much longer than individual 03–05 workflows |
| 01/02 | No reliable benchmark captured in the uploaded final audit | Full runs may require many hours or days; archived intermediate outputs included |
| 08/09 | Typically short formatting/reporting jobs | Operate on archived outputs; do not rerun the underlying catalogue/CBPP calibrations |

The 03/04/05 full-scale run consumed approximately 47.5 minutes of wall time for the slowest successfully completed process (05), aside from recovery/checking. Elapsed time depends on machine, RAM, operating system and R/BLAS/compilation configuration.

## Hardware provenance

Captured on the Mac that generated the verified results: MacBookPro16,1; **32 GiB RAM**; **macOS 26.7.1**; **R 4.6.0**, x86_64 Apple Darwin. See `COMPUTATIONAL_HARDWARE.txt` (hostname removed), `R_SESSION_INFO.txt` and execution logs. Timing values above are empirical observations rather than universal performance guarantees.

## Parallelization and adaptation

- The full-scale empirical controller started **three** separate R processes for scripts 03/04/05, with outputs isolated by workflow. Their **B=1000 simulation loops were serial within each workflow**. Running six independent worker processes would require a documented work split and RNG redesign; simply duplicating workers is not equivalent.
- Script 06 runs Part B with `parallel::mclapply(..., mc.cores=N_CORES, mc.preschedule=FALSE)` when `N_CORES > 1`, otherwise sequential `lapply`. Its source currently sets `N_CORES <- 2`; each replicate uses deterministic seeds derived from scenario and replicate indices, and archived checkpoint `.rds` files permit inspection. On Windows, where `mclapply` multicore forking is unavailable, use `N_CORES <- 1` (serial).
- Script 02's `MAKEFLAGS` setting controls parallel package **compilation**, not statistical simulation parallelism. It detects local CPU count. It should not be mistaken for a 6-worker Monte Carlo run.
- For lower-memory computers, run each workflow sequentially from the root, and use reduced-run modes or archived outputs to avoid redoing long simulations. Avoid changing RNG seeds or statistical settings if comparing against archived reference values.

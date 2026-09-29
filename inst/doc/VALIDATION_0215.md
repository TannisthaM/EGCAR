# EGCAR 0.2.15 validation and limits

Release source: 0.2.14 GitHub commit 3eb1b62e6922868ce2ba7e9916ae1d777b2f5b61,
with the accompanying 0.2.15 patch. Reviewed 28-29 September 2026.
No GitHub push was performed.

## Changes and numerical contract

- EGCAR fit_time now includes optimization through loading extraction.
  solver_time and loading_time are components, not additional charges.
  EGCAR CV total_time = tuning_time + fit_time. Public wall_time includes
  preparation performed in that call. Existing results cannot be corrected
  without rerunning missing timing stages.
- Native group norms and coefficient-error reductions stream edge blocks.
- Direct paths reuse one bounded support-factor cache. Scoring reuses projections.
- Population truth uses the latent-span eigenproblem, retaining marginal
  eigenvalue floors, covariance model and the leading truth subspace.
- EGCAR/SGCA CV now requires converged, finite-scoring folds by default.
  cv_require_convergence=FALSE explicitly restores the older eligibility policy.
- Successful-fit means/plots exclude unsuccessful statuses. Raw failures remain;
  completion_counts.csv counts recorded rows, not missing intended jobs.
- Tier launchers require 0.2.15 and apply the same successful-fit aggregation.

Objectives, penalties, ADMM tolerances, loading rank and SGCA algorithms were not
changed. This release is not an equal-accuracy calibration or an optimality claim.

## Checks completed

- Source builds and installs using R 4.3.3 on Linux; native compilation succeeds.
- R CMD check --no-manual --no-tests --no-vignettes: 0 errors, 0 warnings,
  4 notes: unavailable suggested packages; retained C++14 declaration; installed
  shared-library size; historical LaTeX source under inst/doc. Tests were run
  separately. Final tier-aggregation edits were checked by a separate fixture.
- Regression suite plus focused rerun: 50 test blocks, 921 passing expectations,
  one failed expectation, one errored test block, and one optional test skipped.
  Both outstanding issues are in SGCA subprocess tests; details below.
- All eight new 0.2.15 test blocks pass: 65 expectations, including streamed
  reductions/input nonmutation, dense-versus-reduced truth for K=2..4 and r=1..3,
  deliberately delayed loading timers, convergence eligibility, direct scores,
  warm paths and dense-versus-streamed experiment diagnostics.
- A complete small study (n=35, dimensions 4/5/6, rank 1, two folds) produced
  status=ok for both EGCAR rate and CV methods and both oracles. Added timing and
  completion-count columns were checked. Comparators were explicitly disabled.
- Tier 1/2/3 aggregate fixtures each correctly reduced one success, one timeout
  and one nonconverged fit to a successful-time mean of 1, while preserving all
  three recorded runs and their completion counts.
- All 190 top-level R functions were counted with the R parser and documented;
  27 exports and 8 registered S3 methods remain unchanged.
- The standalone LaTeX guide compiles with pdflatex; its 33-page rendering was
  inspected. No external images/bibliography files are needed.

## Outstanding environment-dependent checks

The SGCA supervisor comparison test failed one expect_warning(..., NA) assertion
because processx/ps reported that a process stat file could not be read. The
process-descendant termination test errored when processx initialization in a
child failed with ps::ps_handle(): Unknown errorfs_error0NA. Seven warnings were
recorded in the full attempt. These checks were not silently skipped or claimed
as passing. Most supervisor numerical/result checks did execute successfully.
Recheck supervision, descendant cleanup and parallel behavior on the cluster.

The optional RGCCA/PMA comparison suite was not enabled. Those packages,
ggplot2 and RhpcBLASctl were unavailable here. Plot aggregation data were checked,
but ggplot2 output rendering and fresh optional-comparator equivalence were not.
No large Midway run, whole-process peak-RSS study or tighter-tolerance stability
calibration was performed. A clean unrestricted full R CMD check is not claimed.

## Focused component benchmark

Elapsed medians after warm-up; Rprofmem allocation totals, not process peak RSS
and not native heap allocation measurements. Exact inputs/results are compared
before measurement. Group/error cases: three views of 700 variables, five repeats.
Population case: three views of 100 variables, rank two, three repeats.

| Operation | Previous seconds | 0.2.15 seconds | Previous R allocation bytes | New R allocation bytes |
|---|---:|---:|---:|---:|
| Incident group norms | 0.013 | 0.002 | 23,588,064 | 19,496 |
| Coefficient error/zero reductions | 0.176 | 0.006 | 182,297,616 | 2,552 |
| Population construction | 0.186 | 0.022 | 21,158,280 | 8,037,360 |

These local component ratios cannot be extrapolated to full solver runtime,
cluster performance or tolerance-matched accuracy. The original 5q-to-3q+O(p)
L21 state compression is retained; it is not newly introduced in 0.2.15.

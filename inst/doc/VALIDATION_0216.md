# Validation of EGCAR 0.2.16

Validated on 2026-09-29 in local Linux, R 4.3.3, g++ 13.3, system BLAS/LAPACK.
This is numerical implementation validation, not a rerun of the user's cluster
pilot or a reproduction of the published simulation tables.

## Outcomes

- Source build and source-archive installation succeeded.
- `R CMD check --no-manual --no-tests --no-vignettes`: **0 errors, 0 warnings,
  4 notes**. Notes concern unavailable optional suggested packages, C++14
  declaration, shared-library size, and the historical installed LaTeX source.
  Unit tests were run separately, with results below.
- **98 new SGCA assertions passed**, in six test blocks. These cover the
  independent Algorithm 1 transcription for ranks 1 and 2; the 2*eta update;
  metric normalization; row sparsity and deterministic ties; fixed-step versus
  optional early stopping; trace-equality Fantope projection; cached/direct
  initializer equality with separate counters; training-size rho scaling;
  common versus paper validation scores in dense and data representations;
  capped-fit/convergence eligibility; study integration and saved diagnostics.
- The broader suite recorded **1019 passing assertions, one failed assertion,
  one errored block, one skipped optional test block, and eight warnings**.
  Both outstanding issues are in the existing process-supervision tests:
  an unexpected process-stat-file warning under a 0.01-second deadline, and
  `processx` failing to load via `ps::ps_handle()` when testing nested-child
  termination (`Unknown errorfs_error0NA`). These also occurred in 0.2.15.
  No claim of a fully passing suite is made. The supervisor implementation was
  preserved; its end-to-end descendant termination still needs validation in
  the target cluster environment.
- Supervised versus in-process SGCA numerical comparisons, direct versus cached
  initialization, and the fold-worker comparison passed before the extremely
  short-deadline warning assertion in that test block.
- All cluster shell scripts passed `bash -n`. Cluster jobs were not submitted.

## Unmodified default-profile run

A simulated rank-two example used n=60 and three views of five variables, five
folds and **all public SGCA default controls**, including the six-hour supervisor.
The valid default sparsity grid was {5,10,15}. All 15 fold fits and the final
refit completed exactly 15000 TGD updates, with finite full-rank loadings. The
final loading satisfied L' B L=I to the asserted 1e-9 tolerance.

The fold initializer counts ranged from 62 to 111. The final initializer took
101 steps and met the 0.005 change criterion. The selected full refit had a final
TGD change about 5.23e-5, exceeding the diagnostic 1e-6 threshold. Accordingly,
its status was `ok` with `completed=TRUE` and `converged=FALSE`. This is deliberate:
fixed-step benchmark completion is not relabeled numerical convergence.
The test is evidence of algorithm execution, not matched statistical accuracy.
The local run also emitted process-stat-file warnings; it returned successfully.

## Included evidence and reproduction

The release's `validation/` directory contains:

- `tests-full.csv` and `full-tests-final.log`: broader suite counts and errors;
- `default-profile.R`, `default-profile.log`, `default-fold-results.csv`, and
  `default-refit-diagnostics.csv`: the default-profile run;
- `package-check.log`, `session-info.txt` and `install-source-archive.log`.

From the extracted release directory after installation:

```sh
Rscript validation/default-profile.R
Rscript -e 'library(egcar); testthat::test_dir("egcar/tests/testthat", reporter="summary", env=new.env(parent=asNamespace("egcar")))'
```

PMA/RGCCA comparisons and ggplot2 rendering were not rerun because those optional
packages are unavailable locally. No large-p peak-memory or cluster-wall-time
benchmark was run. EGCAR's existing 0.2.15 numerical and timing tests remain in
the suite and pass. Use fresh output folders, configurations and CV-data objects
for cluster runs of this release.

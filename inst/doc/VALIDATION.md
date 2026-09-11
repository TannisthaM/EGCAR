> Historical validation report for version 0.2.0. For the current speed-only
> revision, see SPEED_VALIDATION_021.md in this directory.

# Validation of egcar 0.2.0 and the three matching rewrites

## Executed in the creation environment

* Lexical string/delimiter and trailing-comma checks over 52 R files.
  This is not R's parser and is not a claim that `parse()` or `R CMD check` ran.
* Resolution checks for 89 numerical-engine function names and 27 exported
  functions; each export has an R help alias.
* Exact source-text comparisons for 20 retained comparison/simulation functions
  against the user's uploaded local experiment. All matched, including the
  five family CV routines, common CV engine, validation score, penalized SGCA
  TGD, MultiCCA Gram implementation, population generator and support oracle.
* The standalone and package experiment loops match after removing only the
  function wrapper/indentation. Package execution uses installed namespace
  functions; it does not source the standalone script.
* The Rcpp native source, registration source and platform Makevars files are
  unchanged from the preceding accelerated package.
* Repository scans found no removed tied estimator/result labels or separate
  initializer-package name/dependency. The required four initializer functions
  and their copyright/permission notice are included locally.
* Compiled the standalone C++ core with g++ 14.2.0 and system Armadillo; all
  **120 dense/native comparison cases passed**. These include both penalties,
  singular/full-rank covariances, warm starts, zero/large coefficients,
  adaptive augmentation and differing residual-check schedules.
  The maximum reported matrix row-sum norm difference was
  **8.82678305447193e-11**; maximum residual difference was
  **4.99600361081320e-15**. Stopping iterations, convergence flags and final
  augmentation parameters matched in every tested case.

The exact machine-readable source comparisons and native statistics are in
`source_audit.json`. These numerical tests exercise the retained C++ iteration
core, not the Rcpp marshalling layer or newly bundled SGCA code.

## Not executed here

R and Rscript were not available, and network access from the execution
container did not allow installing them. Consequently, R parsing, package
installation, the Rcpp interface, testthat, `R CMD build`, `R CMD check`, SGCA
initializer equivalence in R, multisession CV, plotting and the full experiment
were **not run**. No measured R speedup or clean package-check result is claimed.

## Included checks to run in R

From the extracted package root, install core/check dependencies and run:

```bash
Rscript inst/examples/00_install_dependencies.R --benchmarks --checks
Rscript tools/check.R .
```

Optional comparison and parallel tests are enabled with `EGCAR_TEST_OPTIONAL=true`.
After installation, from the directory containing the delivered standalone script:

```bash
Rscript verify_EGCAR_rewrites.R run_EGCAR_local.R rewrite_checks 1
Rscript verify_EGCAR_rewrites.R run_EGCAR_local.R rewrite_checks_parallel 2
```

The parity checker compares function bodies, folds, CV selections, operator
estimates and statistical metrics for small standalone/package studies under
both the accelerated R and compiled backends. Wall-clock times are intentionally
not compared. Missing RGCCA/PMA dependencies are reported as skipped rather than
passed. The supplied short-iteration smoke configurations can legitimately
produce nonconvergence diagnostics and are not performance studies.

## Deliberate behavior changes and retained definitions

Both EGCAR CV grids are strictly positive and zero is rejected in caller-supplied
grids. The unused coefficient remains zero, as do both coefficients of Oracle1.
The SGCA initializer grid and every other benchmark tuning grid are unchanged.
Combined/tied EGCAR fitting interfaces are absent; the original Oracle1 consensus
splitting is retained privately with both coefficients hard-coded to zero.
All-candidate EGCAR CV failure returns `no_valid_cv`, with no arbitrary coefficient
outside the grid. The old hidden coefficient-one fallback is not retained.
The prior package's missing stats::toeplitz import is fixed.

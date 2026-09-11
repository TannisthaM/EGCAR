# Runnable examples

00_install_dependencies.R: explicit core/benchmark dependency installation.
01_quickstart.R: four EGCAR fits with positive CV grids.
02_all_methods_cv.R: the full original local study through the installed package.
03_parallel_cv.R: shared fold-level parallelism.
04_backend_validation.R: native/R/reference comparison.
05_backend_timing.R: timings on the user's own R installation.
06_comparison_cv.R: all four original common-loss comparison wrappers.

No fitting or test function silently installs a package.

## 0.2.1 single-configuration check

`06_small_all_methods_check.R` uses n=60, p_total=12, r=1, signal=0.8,
five shared folds and five workers. It attempts all ten methods with small
CV grids and records finite loadings separately from convergence. Use this
execution check before a large study. It was supplied but not run in the
R-less preparation environment.

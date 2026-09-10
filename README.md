# egcar 0.2.1

Accelerated sparse generalized correlation analysis via pairwise regression.
This revision provides **separate L11-only and L21-only estimators**, each with
positive-grid cross-validation and rate-scaled fitting. It also includes the
complete local experiment and its original common-loss comparison CV routines.

## Speed-only revision

The compiled EGCAR solver adds mapped Eigen products for small matrices,
read-only mapped covariance inputs, reusable matrix workspaces and fused
proximal/dual/residual loops. Large products retain Armadillo/BLAS. This extends
the computation ideas in `ZixuanWu1/EfficientCCA` without adopting its stopping
rule, its two-view loss, or its group-penalty definition. See
[the source review](inst/doc/EFFICIENTCCA_REVIEW.md) and
[validation details](inst/doc/SPEED_VALIDATION_021.md).

**Every `R/` implementation file is unchanged from 0.2.0.** In particular,
SGCA, RGCCA, SGCCA, MultiCCA, the oracles, statistical objectives, CV grids,
shared folds and worker allocation are not changed. Use `backend = "cpp"`
for the new native implementation. `backend = "R"` and `"reference"` retain
the prior implementations. Install the new `RcppEigen` build dependency.

A one-configuration check of all ten methods, with five CV workers, is:

```r
source(system.file("examples", "06_small_all_methods_check.R", package = "egcar"))
```

The check uses n=60, total dimension 12, rank 1, signal 0.8 and small positive
EGCAR CV grids. It is an execution diagnostic, not a performance comparison.

## Repository layout

Place the contents of this `egcar/` folder at your GitHub repository root.
`DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, and `man/` should be immediately visible.
The method-specific `R/egcar.r`, `R/reduced_rank_regression.R`,
`R/group_reduced_rank_regression.R`, `R/alt_*.R`, helpers and metrics follow the
organization of https://github.com/cran/ccar3/tree/master/R . The compiled
`src/` directory is intentionally retained for the accelerated native backend.

## Install after publishing through GitHub Desktop

In RStudio, replace `YOUR_GITHUB_USERNAME/egcar` with your actual owner/repository:

```r
install.packages(c("remotes", "Rcpp", "RcppArmadillo", "RcppEigen", "future", "future.apply",
                   "RGCCA", "PMA", "ggplot2", "RSpectra", "RhpcBLASctl"),
                 repos = "https://cloud.r-project.org")
remotes::install_github("YOUR_GITHUB_USERNAME/egcar", dependencies = NA,
                        upgrade = "never", build_vignettes = FALSE)
library(egcar)
packageVersion("egcar")
```

The package name used by `library()` is lowercase `egcar`, regardless of the
GitHub repository's capitalization. Source installation needs an R-compatible
C++14 build toolchain. On Windows use Rtools matching your R version; on macOS
use the appropriate command-line build tools. A separate Armadillo installation
is not needed for the R package: RcppArmadillo supplies its headers.
RcppEigen supplies the Eigen headers for the small-matrix native branch.

The SGCA initializer and its full four-function dependency closure are bundled.
No separate SGCA implementation package is needed or installed. RGCCA supplies
the RGCCA and SGCCA fitters, and PMA supplies MultiCCA and its threshold helpers.
The comparison CV routines are the ones from the user's local experiment, NOT
alternative tuning routines from the underlying libraries or source repository.

## Four EGCAR variants

```r
sim <- egcar_simulate(n = 120, p_list = c(15, 15, 15), rank = 2,
                      active_per_view = 5, seed = 12)
shared <- egcar_cv_data(sim$views, nfolds = 5, seed = 13)
ctl <- egcar_control(backend = "cpp")
l11_cv <- EGCAR_L11_CV(shared, rank = 2, lambda = 10^seq(-5, 4), control = ctl)
l21_cv <- EGCAR_L21_CV(shared, rank = 2, lambda = 10^seq(-5, 4), control = ctl)
l11_rate <- EGCAR_L11_Rate(shared, rank = 2, control = ctl)
l21_rate <- EGCAR_L21_Rate(shared, rank = 2, control = ctl)
```

`EGCAR_L11()` and `EGCAR_L21()` fit at a specified `lambda`. The generic
`egcar_fit`, `egcar_cv`, `egcar_rate`, and `egcar_path` interfaces remain available.
L11 adds only the entrywise penalty; L21 adds only the sum of concatenated-row
Euclidean norms across all incident edges. L21 iterates C/G/V without an entrywise
proximal step, and its linear solve retains the factor `2 * mu`.

**Both CV functions reject zero.** Their defaults are `10^seq(-5, 4)`.
The unused penalty is zero by construction, not an unregularized CV candidate.
Zero remains allowed for fixed fits and the original zero-penalty population
oracle. SGCA's separate initializer grid is unchanged and can include zero.

The rate scales are `sqrt(log(p)/n)` for L11 and
`sqrt((d_max + log(p))/n)` for L21, as in the supplied experiments.
The second is retained as a benchmark scaling; this package makes no new claim
that it is an optimal statistical rate.

## Run the SAME full local experiment

```r
cfg <- egcar_experiment_config()
ans <- run_egcar_experiments("egcar_package_outputs", workers = 5, n_reps = 1,
                             config = cfg, backend = "cpp")
```

Defaults: three views of size 15; n = 30,45,100,1000,5000,10000;
ranks 1,2,5; five active variables per view; Toeplitz parameters 0.5,0.7,0.9;
signal 0.8; master seed 20260907; five common folds; nested observations as n grows.
The complete local simulation loop and plotting routines are shared with
`run_EGCAR_local.R`, supplied separately and under `inst/standalone/`.
`inst/examples/02_all_methods_cv.R` runs the installed-package version.

All ten benchmark labels:
`Oracle1-population`, `Oracle2-support`, `EGCAR-L11-rate`, `EGCAR-L11-CV`,
`EGCAR-L21-rate`, `EGCAR-L21-CV`, `SGCA`, `RGCCA`, `SGCCA`, `MultiCCA`.
There is no combined-penalty or tied EGCAR method.

Each method gets the same fold-level worker budget, capped by the fold count.
Families run sequentially; no nested method/fold parallelism is introduced.
The numerical kernels retain cached spectral bases and denominators, thin SVDs
for wide views, dimension-dependent multiplication order, full null-space
corrections, warm paths, optional checked partial loading eigensolves, and C++
updates. The native core is compiled on package installation, not on every fit.
No wall-clock speedup in R is claimed without running the timing example locally.

The original Oracle1 zero-coefficient consensus splitting is retained privately
for baseline reproducibility. It is not a public combined-penalty estimator.
Oracle2 still uses sample covariance restricted to the true supports.

## Comparison CV

`sgca_cv`, `rgcca_cv`, `sgcca_cv`, and `multicca_cv` accept the same shared
`egcar_cv_data` object. They preserve the local experiment's common generalized
Rayleigh loss, training-mean centering, covariance divisor n, tie rules, training-only
sign synchronization, cached initializer, penalized TGD and PMA Gram backend.
See `inst/examples/06_comparison_cv.R` for explicit calls and grids.

## Output

Use a new output directory for each run: checkpoints overwrite files of the
same names. Main checkpoint: `egcar_simulation.rds`.
Also saved: `simulation_results.csv`, `cv_grid_results.csv`,
`selected_cv_parameters.csv`, `egcar_cv_fold_results.csv`,
`external_cv_fold_results.csv`, worker/backend metadata, failures, session info,
metric plots with and without Oracle1, and loading-comparison data/PDFs.
Loading graphics use the original post-fit global normalization and Procrustes
alignment to truth; no truth is used in CV or ordinary fitting.
The new fit keys are `EGCAR_L21_rate` and `EGCAR_L21_CV`; L11 keys are unchanged.

Failed CV is recorded as `no_valid_cv`, with no selected penalty or arbitrary
fallback fit. Finite nonconverged fits remain eligible exactly as in the local
benchmark, with convergence diagnostics retained. Smoke tests intentionally use
short limits and can report nonconvergence; they are not accuracy/runtime studies.

## Check before the full study

```r
run_egcar_experiments("egcar_smoke", workers = 2, smoke_test = TRUE)
```

From the package root, with testthat installed:

```bash
Rscript tools/check.R .
```

Set `EGCAR_TEST_OPTIONAL=true` to include optional comparison/parallel tests.
The supplied standalone/package equivalence checker is in
`tools/verify_EGCAR_rewrites.R` and is also provided separately.

**Validation limitations:** R was unavailable in the creation environment.
The Rcpp interface, installation, R CMD check and experiment runs were not
executed there. Static checks and 120 standalone native numerical comparisons
passed. See `inst/doc/VALIDATION.md` for exactly what was and was not tested.

## Attribution and publication metadata

Bundled initializer functions were obtained from
https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R .
The MIT copyright/permission notice for Claire Donnat is retained in the source
and `inst/COPYRIGHTS`. Required functions are `Soft`, `updatePi`, `updateH`, and
`sgca_init_fixed`; unrelated upstream CV/TGD routines are not imported.
Acceleration ideas were checked against ccar3's `R/ecca.r`. The EGCAR native
solver is the existing independently written implementation.

The overall package retains its GPL-3 license; bundled components retain their
notices. Author/maintainer placeholders in DESCRIPTION must be replaced before
public distribution. This archive is an editable source repository, not a
Windows binary package or a claimed CRAN release.

## Standalone versus package validation

```bash
Rscript verify_EGCAR_rewrites.R run_EGCAR_local.R rewrite_checks 1
Rscript verify_EGCAR_rewrites.R run_EGCAR_local.R rewrite_checks_parallel 2
```

This verifier checks source-body identity for the comparison CV routines and
runs small local/package studies with both accelerated R and native backends.
It compares folds, selected coefficients, operator estimates and statistical
metrics, but deliberately does not compare wall-clock times. It reports missing
comparison dependencies as skipped, not passed. It has not been executed in the
creation environment because that environment did not contain R.

# EGCAR 0.2.16: paper-based SGCA

SGCA now implements Gao and Ma's Algorithm 1. Read
[the SGCA specification](inst/doc/SGCA_PAPER_0216.md) before comparing with older results.
Defaults: fixed 15000 updates, eta=0.001, lambda=0.01, no added metric ridge;
rho=0.5*sqrt(log(total p)/training n) is recomputed per fold and refit.
The initializer enforces trace=rank and uses tolerance 0.005 with a real cap of
1000 outer iterations. The initializer tolerance/cap come from authors' R code,
not a paper accuracy guarantee.

Five-fold CV selects total row sparsity on valid members of 5,10,...,100 using
the paper's covariance-trace score. `benchmark_control(sgca_cv_score="common_loss")`
selects the common comparison score. `rho_grid=NULL` is the rate rule; explicitly
supplied numeric rho values remain direct coefficients.

`status="ok"` reports completion with valid loadings. `converged` is a separate
numerical diagnostic: completing 15000 steps does not prove convergence. See
`fit$diagnostics` and `fit$cv_fold_table`. A capped initializer remains usable in
the default fixed-iteration benchmark; `sgca_cv_require_convergence=TRUE` makes
both numerical criteria mandatory. The six-hour complete-CV budget remains.
EGCAR's full-loading runtime accounting and memory improvements from 0.2.15 remain.

```r
library(egcar)
shared <- egcar_cv_data(views, nfolds=5, seed=1)
fit <- sgca_cv(shared, rank=1)
fit$status
fit$converged
fit$diagnostics
```

Reinstall 0.2.16, restart R and workers, copy the updated `inst/cluster` launchers,
and use a fresh result directory. Old launchers explicitly override SGCA grids,
ridge and caps. Versioned documents below describe historical releases; this
section and `SGCA_PAPER_0216.md` govern current SGCA behavior.

# EGCAR 0.2.15 update

This release records the full EGCAR fit through loading extraction. `fit_time`
includes `solver_time` and `loading_time`; never add `loading_time` again.
For CV, `total_time = tuning_time + fit_time`. `wall_time` measures the current
public fit/CV call, including preparation performed in that call. Preparation
stored in a reusable object is historical and must not be charged twice.
Experiment results also report shared preparation separately.

New efficiencies: streaming native group norms and coefficient-error reductions,
shared bounded loading-factor caches across a regularization path, direct reuse
of prediction scores in `egcar_score`, and a reduced population-truth eigenproblem.
Objectives, statistical penalty scales, residual tolerances, and requested loading
rank are unchanged. In 0.2.15 CV required convergence for both methods; SGCA now has the separate benchmark policy described above.
Successful-fit plots exclude unsuccessful statuses; `completion_counts.csv`
reports success/timeout counts among recorded runs. Inspect that denominator
alongside plots. These changes do not establish equal accuracy across methods.

See `inst/doc/VALIDATION_0215.md` for measured savings, tests and remaining validation limits.
Use a fresh study output directory; updated tier launchers require 0.2.15.

Reinstall from source and restart R/workers: this release adds native entry points.
Older saved objects retain their old timing interpretation; use `timing_version`
to distinguish newly computed results.

# Calibrated L21 rate in 0.2.14

`EGCAR_L21_Rate()` and `egcar_rate(penalty="l21")` now use **`multiplier=0.05`**,
selected on independent pilot simulations and frozen before confirmation.
`egcar_experiment_config()$rate_c_g` uses that same multiplier. The rate formula,
objective, solver tolerances, rank and direct-penalty CV are preserved.

```r
cfg <- egcar_experiment_config()
cfg$rate_c_g  # 0.05; also recorded in each experiment CSV
# fit <- EGCAR_L21_Rate(views, rank = 2)
# fit$rate_multiplier
# fit$rate_zero_certified
```

**Remove `rate_c_g = 1` from old experiment scripts to use the new default.**
An explicit `rate_c_g=1` or `multiplier=1` deliberately retains the old rule.
No automatic penalty reduction or rank substitution is performed. A zero or
rank-deficient fit is recorded as `status="invalid_loading"` and
`loading_valid=FALSE`, while `converged` continues to describe the optimizer.
Its unavailable subspace errors remain NA. Zero-solution certificates and
selected-row counts make these cases visible in result and failure CSVs.

The complete pilot protocol, candidate scores, independent confirmation and
limitations are in `inst/doc/L21_RATE_CALIBRATION_0214.md`, with raw tables in
`inst/validation/l21_rate_0214/`. A valid loading alone does not establish good
subspace recovery or sparsity. Calibration is specific to the supplied model
and tested regimes; the coefficient remains configurable.
The selected value is the smallest tested candidate. It produced valid,
converged loadings in all 60 independent confirmation cases (old value: 9/60),
but selected every variable. High-dimensional subspace recovery remained poor.
It repairs over-penalization, but does not establish sparse recovery or an
optimal coefficient outside the tested grid.

Updated tier-1/2/3 R and SLURM launchers are in `inst/cluster/`. They require
0.2.14 and inherit its rate constant. They also use the existing six-hour SGCA
policy instead of overriding it with the old finite iteration caps. Their
signal remains 0.8, as in the supplied scripts; set it deliberately for other
studies. Copy these launchers into your cluster experiment directory and use
a fresh output directory, so old `.done` markers do not skip the rerun.

To install the source ZIP, restart R and run:

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "RSpectra", "callr"))
unzip("egcar_0_2_14_l21_calibrated.zip", exdir = "egcar_0214")
install.packages("egcar_0214/egcar", repos = NULL, type = "source")
packageVersion("egcar")  # 0.2.14
```

# SGCA time budget in 0.2.13

```r
benchmarks <- benchmark_control(sgca_time_limit = 6 * 60 * 60)
fit <- sgca_cv(shared, rank = 1, benchmarks = benchmarks)
fit$status
fit$converged

cfg <- egcar_experiment_config(sgca_time_limit = 6 * 60 * 60)
```

The six-hour budget covers **one whole SGCA CV grid plus final refit**, including
initialization. It is not six hours per candidate or fold. Shared data/fold
preparation happens before the budget starts. In the experiment runner, each
dataset/rank/replicate receives a fresh budget. In 0.2.13--0.2.15 both caps defaulted to `Inf`. Version 0.2.16 uses the fixed-step policy above.
A timeout returns `status="time_limit"`, `timed_out=TRUE`, `converged=FALSE`, and
`L=NULL`, with the reason "SGCA did not converge in real time". This means the
requested procedure did not finish within the budget, not that the optimizer
could never converge. Incomplete fits and CV results are not scored as successes.

A supervised R subprocess can be stopped during a long compiled matrix operation.
Startup/input transfer/cleanup add overhead; the subprocess also holds its own
input copy, so the timeout feature is not a memory optimization. The controller
keeps the requested fold-worker count. Set `sgca_time_limit=Inf` to run in-process
without supervision. For finite iteration limits, exhaustion remains
`status="not_converged"`, not `"time_limit"`.

See `inst/doc/egcar_efficiency_writeup.tex` for the full derivation and audit.

# egcar 0.2.12

Accelerated sparse generalized correlation analysis via pairwise regression.
This revision provides **separate L11-only and L21-only estimators**, each with
positive-grid cross-validation and rate-scaled fitting. It also includes the
complete local experiment and its original common-loss comparison CV routines.


## Algorithm implementation (0.2.12)

The default L21 solver stores `C`, `Hk`, `Hl`, and per-row multipliers `a`.
It reconstructs the group and scaled-dual values only while processing each
edge. Main ADMM state is `3q + p` rather than `5q`, where
`q = sum(p_k * p_l)` over unordered view pairs. This is a state-memory reduction,
not a claim of 40% lower total process RAM.

Both EGCAR penalties recover the leading requested eigenpairs through a
matrix-vector callback over the edge blocks. Loading normalization uses compact
spectral factors; the selected symmetric operator and its dense square roots
are not assembled. `RSpectra` is now required. Failed partial solves produce an
explicit error after a larger-subspace retry, without a hidden dense fallback.
A two-variable operator uses an exact analytic leading pair.

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "RSpectra"))
install.packages("egcar_0.2.12.tar.gz", repos = NULL, type = "source")
library(egcar)
control <- egcar_control()  # compact state and partial eigenpairs enabled
# fit <- egcar_fit(views, rank = 2, penalty = "l21", lambda = 0.02, control = control)
```

Old fits and G/V or endpoint warm-start lists remain accepted. For code that
needs explicit old-style `fit$solver$G` and `$V`, use
`egcar_control(compact_state = FALSE)`. This affects retained output only.
Restart R and any CV workers after installation because the native state API
has changed. `keep_full_C=TRUE` and `partial_eigen=FALSE` explicitly request
dense diagnostic output/computation; the reference backend is also dense.

See [ALGORITHM_IMPLEMENTATION_0212.md](ALGORITHM_IMPLEMENTATION_0212.md) for
implementation details, numerical checks, timings and the remaining memory costs.

## Speed and memory update (0.2.11)

This release removes redundant dense preparation, reuses full-support loading
spectra, computes wide-fold validation scores through projected observations,
and reduces L21 working arrays and CV state conversions. Existing public calls, penalties, default tolerances, grids and centering rules
are retained. The new optional control `loading_cache_max_bytes` defaults to 64 MiB.
The 0.2.10 future-plan fix remains in place.

See [OPTIMIZATION_REPORT.md](OPTIMIZATION_REPORT.md) for the changes, ccar3 review,
exact memory calculations, validation limitations and benchmark instructions.
At the time of the 0.2.11 release, only algebra and extracted-kernel checks
were available. The 0.2.12 checks use an installed R package; consult its report
for current validation results.

After installing the extracted source with `R CMD INSTALL /path/to/egcar`, run:

```r
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.12")
source(system.file("examples", "08_performance_check.R", package = "egcar"))
```

The historical `inst/standalone/run_EGCAR_local.R` is retained as a comparison
snapshot; use the installed package for these optimizations.

## Memory-first experiment storage (0.2.3)

`run_egcar_experiments()` now defaults to compact storage: it does not retain full solver/package fits, fold-level CV diagnostics, or full loading-plot data. It saves the main result/CV-summary tables and small xz-compressed per-configuration loading archives under `compact_loadings/`. Restore diagnostic payloads with `egcar_experiment_config(save_fits = TRUE, save_cv_fold_results = TRUE, save_loading_data = TRUE, retain_benchmark_fits = TRUE)`. These switches affect retention only; the estimators, grids, folds, seeds, scoring rules, and parameter selection are unchanged.


## Matrix-interface patch (0.2.2)

The native return path now allocates each state block as an explicit R numeric
matrix instead of relying on generic nested-container wrapping. The R dispatcher
checks the native API marker and every edge shape before row-norm calculations.
This targets the reported `rowSums()` dimension error; the exact original failing
call was not reproduced in the patching environment, where R was unavailable.
The numerical solver core and experiment settings are unchanged.

After a clean source reinstall and a fresh R session, run the installed check:

```r
source(system.file("examples", "07_matrix_interface_check.R", package = "egcar",
                   mustWork = TRUE))
```

It checks native matrix dimensions/values, requested ranks 1, 2 and 5 at signal
0.8, both separate penalty families, warm starts, paths, rates, and serial versus
five-worker CV. This is a diagnostic with small grids, not the full experiment.
Set `EGCAR_MATRIX_CHECK_WORKERS=1` for a serial-only check. Details and exact
validation limitations are in [the patch notes](inst/doc/MATRIX_INTERFACE_FIX_022.md).

## Previous speed-only revision (0.2.1)

The compiled EGCAR solver adds mapped Eigen products for small matrices,
read-only mapped covariance inputs, reusable matrix workspaces and fused
proximal/dual/residual loops. Large products retain Armadillo/BLAS. This extends
the computation ideas in `ZixuanWu1/EfficientCCA` without adopting its stopping
rule, its two-view loss, or its group-penalty definition. See
[the source review](inst/doc/EFFICIENTCCA_REVIEW.md) and
[validation details](inst/doc/SPEED_VALIDATION_021.md).

**In revision 0.2.1, every `R/` implementation file was unchanged from 0.2.0.** In particular,
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
`egcar_cv_data` object. SGCA defaults to the paper covariance-trace CV score;
the other wrappers retain the common generalized Rayleigh loss. All use the
same fold assignment and training-mean centering. EGCAR 0.2.16 changes SGCA's
solver, tuning defaults and eligibility as described at the top of this file.
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
fallback fit. EGCAR requires convergence and finite requested-rank scores on every fold. SGCA uses its own completion/convergence controls described above. Smoke tests intentionally use
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
The former standalone-equivalence checker now runs the maintained numerical
conformance tests: `Rscript tools/verify_EGCAR_rewrites.R .` from the package root.
The standalone entrypoint delegates to the same installed package.

**Validation history:** Early releases were checked without an R runtime.
Version 0.2.12 was subsequently compiled and tested in R: 767 expectations
passed and `R CMD check --no-manual` had zero errors/warnings. Its performance
measurements use the dedicated benchmark script, not the full all-methods
example or SLURM sweep. See `ALGORITHM_IMPLEMENTATION_0212.md` and the new
implementation write-up for scope and version-specific validation.

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

## Numerical conformance checks

From the extracted package source directory, after installing 0.2.16:

```sh
Rscript tools/verify_EGCAR_rewrites.R .
```

This runs the independent SGCA Algorithm 1 checks plus EGCAR numerical and timing
regressions. The standalone entrypoint now delegates to the installed package;
there is no separate source-body identity claim. See `inst/doc/VALIDATION_0216.md`
for actual validation outcomes, including the local process-supervision limits.

# EGCAR 0.2.11: time and memory changes

Updated from the supplied `egcar_0_2_10_plan_nse_fix.zip` (version 0.2.10).
The implementation targets both EGCAR penalties, fixed/rate fits, cross-validation,
and the installed-package experiment runner. No additional R dependencies are
introduced. Existing public calls remain valid. The 0.2.10 future-plan/NSE fix,
short cache keys, dimension checks and overflow fixes are retained.

**Validation status:** independent numerical algebra checks and extracted C++
group-kernel checks passed. R, Rcpp and Armadillo were unavailable in the editing
environment, and dependency downloads were blocked. Consequently, installation,
R parsing, `R CMD check`, the complete native solver, parallel R execution and
end-to-end R timings have **not** been verified. No measured package speedup or
peak-RAM reduction is claimed. Run the included R checks before a large study.

## Changes

| Area | Supplied package | Version 0.2.11 | Effect |
|---|---|---|---|
| Wide training views | Full covariance eigendecomposition, then a separate thin data SVD | Thin SVD directly; same positive-singular-value rule | Removes an unused cubic-cost factorization and large eigenvector array |
| Prepared objects | Always stores total-dimension block-diagonal covariance and full edge eigenvalue products | Builds these only for dense reference/oracle solves | Removes redundant quadratic caches from ordinary EGCAR preparations |
| Validation when `n_validation < p_total` | Stores both dense total covariance and dense block-diagonal covariance | Stores centered validation views and contracts through component scores | Uses `n_validation * p_total` instead of `2 * p_total^2` doubles |
| Loading normalization | Repeats covariance eigenanalysis even when every row of a view is selected | Reuses that view's cached spectrum; retains the ridge on thin-basis complements | Removes repeated full-support eigenanalysis |
| Factorization cache | At most four entries, irrespective of their size | Both entry-count and 64 MiB default byte limits | Prevents retained factorization caches from growing to multiple GB |
| Native L21 iterations | Retains `Wk = C + Vk` and `Wl = C + Vl` for all edges | Accumulates norms directly and recomputes each scalar during the proximal update | Removes two full edge collections and their write/read traffic |
| Transformed coefficients | Retains a matrix per edge in every solve | Uses one workspace except when group objective history needs all edges | Lowers working storage |
| L21 CV warm starts | Assembles view-wide G/V matrices and slices them back into endpoint matrices at the next coefficient | Carries endpoint matrices directly between candidates | Avoids repeated large reshaping/copying |
| EGCAR parallel folds | Task closures can retain the entire calling frame, including all folds | Small explicit closure environments; unused training views omitted from EGCAR payloads | Reduces avoidable worker serialization and retention |

The R backend also avoids persistent Wk/Wl edge collections. Public fitted
objects still expose their established G/V warm-start fields. Both legacy and
compact internal warm starts are accepted. Shared CV-data objects continue to
retain training views for SGCA/RGCCA/SGCCA/MultiCCA comparisons.

The comparison estimators, oracle definitions, covariance divisor, penalty
coefficients, grids, training-fold centering, ADMM residual criteria, iteration
limits, adaptive-mu rules and positive-loading-rank requirements are unchanged.
Floating-point evaluation order changes, so equivalence is numerical, not
bit-for-bit. Zero-ridge and partial-support loading factorizations retain the
original dense calculation.

Two existing validation inconsistencies were also addressed: a missing native
API marker now triggers the already-documented error, and the core centering
check no longer accesses fold indices/means removed by an earlier release.

## Why the score and loading changes preserve the formulas

For training-centered validation views X_k and loading blocks L_k, define
T_k = X_k L_k. The two score matrices are exactly

- A = crossprod(sum_k T_k) / n_validation;
- Q = sum_k crossprod(T_k) / n_validation.

These equal `L' Sigma L` and `L' Sigma0 L`. The existing ridge, eigenvalue floor
and normalized trace calculation are then applied without modification.

For a full-support covariance with thin spectrum S = Q diag(d) Q', the regularized
root is `sqrt(a) I + Q diag(sqrt(d+a)-sqrt(a)) Q'`, where `a` is the existing ridge.
The implementation applies the existing eigenvalue floor to both the retained
spectrum and its complement, and uses the corresponding inverse-root identity.
It does not discard the nullspace ridge. The ADMM solve likewise retains its
original full nullspace contribution.

## Analytical memory example, not a measured peak

Consider three views of 1,000 variables, 200 observations and five balanced
folds (40 validation rows each). Counting only matrix payloads:

| Removed or reduced allocation | Bytes saved | MiB saved |
|---|---:|---:|
| Dense reference-only caches in full + five fold preparations | 576,000,000 | 549.3 |
| Five pairs of dense validation covariances, replaced with validation views | 715,200,000 | 682.1 |
| Two L21 work arrays, per active native solve | 48,000,000 | 45.8 |

The first two rows sum to approximately **1.20 GiB of retained matrix payload**,
excluding savings from the discarded wide-view covariance eigenvectors. These
counts omit R object overhead, temporary allocations, sharing, allocator reuse
and worker behavior. They are not process-RSS measurements.

Dense marginal/cross-view covariance blocks, coefficient matrices and localized
loading matrices still require substantial memory. This release does not make
the estimator fully matrix-free. In particular, very large selected supports
still require dense loading products/eigensolves, and running more workers
still multiplies the mutable solver working set.

## ccar3 review

Sources inspected on 2026-09-22:

1. [CRAN ccar3 manual, version 0.1.2](https://cran.r-project.org/web/packages/ccar3/ccar3.pdf).
2. [Distributed R/ecca.r source](https://rdrr.io/cran/ccar3/src/R/ecca.r).
3. [Distributed R/helpers.r source](https://rdrr.io/cran/ccar3/src/R/helpers.r).
4. [Distributed R/reduced_rank_regression.R source](https://rdrr.io/cran/ccar3/src/R/reduced_rank_regression.R).

The available source mirror identifies itself as ccar3 0.1.0, older than the
current manual; the current CRAN source archive was not retrievable here.

`ecca_across_lambdas` reuses spectral quantities across penalties. EGCAR already
used that idea; this revision removes additional preparation and warm-start
conversion costs. The source also scores through projected observations, which
motivated the multiview score contraction above. Its helpers use mapped compiled
multiplication, already present in EGCAR. The RRR code localizes postprocessing
to active rows; EGCAR already does this and now reuses full-support spectra.

The newer manual describes matrix-free conjugate-gradient solves for RRR.
Those are not transplanted into EGCAR's different bilinear covariance operator;
its existing exact spectral solve is retained. Likewise, ccar3's CV covariance
downdates/global preprocessing and two-view scoring are not substituted for
EGCAR's training-fold centering and multiview criterion. This is an independent
implementation; no ccar3 source code was copied.

## Checks and installation

After extracting the zip, install the `egcar` directory in a fresh R session:

```bash
R CMD INSTALL /absolute/path/to/egcar
Rscript /absolute/path/to/egcar/tools/check.R /absolute/path/to/egcar
```

`tools/check.R` requires the existing test/build dependencies, including testthat,
and runs R parsing, build and package checks. Optional comparison-package tests
remain optional. A shorter installed-package diagnostic is:

```r
library(egcar)
stopifnot(packageVersion("egcar") == "0.2.11")
source(system.file("examples", "08_performance_check.R", package = "egcar"))
```

To benchmark against the original source, keep its extracted package directory
separate and run:

```bash
Rscript /path/to/revised/egcar/tools/benchmark_0211.R \
  /path/to/original/egcar /path/to/revised/egcar /path/to/benchmark_results 3
```

The driver installs versions in separate local libraries, runs separate R
processes, alternates version order, and uses identical data, folds, penalties,
thread counts and fixed iteration budgets. It writes timings, object/serialized
sizes, numerical equivalence checks and session information. Its default cases
are diagnostics; edit the worker's case list to match the study dimensions.
Object and serialized sizes are explicitly distinguished from peak process RAM.

Dependency-light checks actually executed here:

| Check | Cases | Maximum relative discrepancy |
|---|---:|---:|
| Dense versus data-space score, independent NumPy formulas | 48 | 2.70e-9 |
| Dense versus cached loading-root formulas | 32 | 1.27e-11 |
| Full nullspace ADMM linear equation | 12 | 5.64e-15 |
| Extracted old/new C++ group kernels, including residuals | 56 | 8.89e-16 |

The C++ kernel check compiles the changed functions verbatim with a small owning
matrix facade; it does **not** compile Rcpp/Armadillo or the complete native
solver. Additional source checks cover balanced delimiters, native-header
synchronization, exported names and the retained future-plan fix. These are not
an R parser or substitute for `R CMD check`.

Reproduce the algebra/kernel checks with Python + NumPy/SciPy and a C++14 compiler:

```bash
python tools/check_algebra_0211.py
python tools/check_group_kernels_0211.py
```

Results from this editing session are under `inst/validation/performance_0211_*`.
The full R regression tests are in `tests/testthat/test-performance-0211.R` and
cover wide/mixed/single-column views, dense-reference agreement, deficient score
matrices, warm starts, cache limits, CV selection and multisession payloads.

The independent historical `inst/standalone/run_EGCAR_local.R` is deliberately
retained as a comparison snapshot. Use `egcar::run_egcar_experiments()` for the
optimized experiment runner. `tools/verify_EGCAR_rewrites.R` compares its score
numerically and continues to compare unchanged comparison-method code.

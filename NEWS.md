# egcar 0.2.16

- Replace the SGCA backtracking variant by paper Algorithm 1 with fixed steps,
  correct gradient factor, initial/final metric normalization and row thresholding.
- Use trace-equality Fantope projection and leading algebraic eigenvectors.
- Set fixed 15000-step TGD, lambda 0.01, eta 0.001, zero metric ridge, and
  fold-specific rho=0.5*sqrt(log(p)/n_train). Enforce the initializer cap of 1000.
- Default SGCA CV to the paper score and sparsity grid; keep common scoring optional.
- Separate benchmark completion from numerical convergence; retain compact
  per-stage diagnostics and actual rho in results and fold tables.
- Update cluster launchers so earlier explicit settings do not override these defaults.
- Preserve 0.2.15 EGCAR timing through loading extraction and memory improvements.

# egcar 0.2.15

* Include loading extraction in public and experiment EGCAR fit times; retain
  solver/loading components and avoid double counting in CV totals.
* Stream group norms and operator errors from edge blocks in native code.
* Reuse bounded support factor caches across direct regularization paths.
* Reuse projected scores for scoring; reduce simulation truth eigenproblems.
* Require converged EGCAR/SGCA CV folds by default (explicit legacy opt-out).
* Exclude failed rows from successful-fit summaries and save completion counts.
* Update tier launchers to require 0.2.15 and filter unsuccessful plot rows.
* Add numerical equivalence, timer-boundary, and eligibility regression checks.

# egcar 0.2.14

- Calibrate one fixed L21 rate multiplier using 36 independent pilot scenarios
  and held-out scores. Freeze it before a separate 60-scenario confirmation.
  Ship the protocol, scripts, raw tables and limitations. No supplied experiment
  errors or population truths are used to select the multiplier.
- Use the frozen coefficient 0.05 in the default L21 rate API and experiment config.
  Explicit numeric multipliers keep their original meaning; L11 and direct CV
  grids retain their defaults. The objective and rate formula are unchanged.
- The selected coefficient is the smallest tested candidate and produces dense
  support. Independent confirmation returned 60/60 valid, converged loadings
  compared with 9/60 for the old coefficient. High-dimensional subspace recovery
  remained poor; the report describes this limitation explicitly.
- Distinguish optimizer convergence from loading validity in experiment CSVs.
  Zero or deficient-rank rate fits now report invalid_loading with a reason.
  Preserve external error/timeout/skip statuses and keep undefined metrics NA.
- Add a tiled, covariance-block zero-solution certificate and selected-row,
  zero-estimate and loading-validity diagnostics. No full matrix is constructed
  for the certificate, and it never silently changes a penalty or rank.
- Include updated tier launchers requiring 0.2.14, inheriting the calibrated
  rate and using the six-hour SGCA policy from 0.2.13.

# egcar 0.2.13

- SGCA now uses a six-hour wall-clock budget for one complete common-loss CV
  grid and full-sample refit, including initialization. The budget is shared by
  folds/candidates and is reset for each dataset/rank/replicate.
- `benchmark_control(sgca_time_limit=21600)` and
  `egcar_experiment_config(sgca_time_limit=21600)` configure the budget.
  Both SGCA iteration caps default to `Inf`; explicit finite caps still work.
- A supervised `callr` process terminates the SGCA process tree on timeout,
  including compiled matrix operations. Timeout returns `status="time_limit"`,
  `converged=FALSE`, and no loadings. Shared fold/data preparation is outside
  the budget. Startup, serialization and cleanup add overhead and another
  process can increase peak memory.
- Added a detailed mathematical implementation note with LaTeX source. It
  separates implemented EGCAR changes from proposed comparison-method changes,
  and states exactly which benchmark scripts were run for version 0.2.12.

# egcar 0.2.12

* Implement the compressed L21 ADMM representation C/Hk/Hl/a in both native
  and R solvers. Group norms are accumulated during the C pass, followed by
  one proximal/dual pass. Main dense state decreases from 5q to 3q+p scalars.
* Preserve full-space residuals, objective histories, null-space corrections,
  adaptive-mu dual rescaling, and arbitrary legacy warm starts. Legacy state
  becomes compressible after its first ordinary proximal update.
* Return compressed warm-start state by default. Set compact_state=FALSE for
  legacy view-wide G/V output. Cross-validation always keeps compact state.
* Recover only the leading r algebraic eigenpairs using an RSpectra function
  callback and a native selected-edge matrix-vector product. No full selected
  operator or dense covariance square-root matrices are constructed by default.
* Store loading metric factors as bases and scalar spectral weights, retaining
  the established ridge scaling, metric floor, and entire covariance range.
* Require RSpectra. Partial solves retry with a larger Krylov space and never
  silently fall back to a dense EVD. A 2x2 problem uses its analytic leading
  pair. Explicit partial_eigen=FALSE and rank equal to dimension permit a full
  solve; the dense reference backend remains available for diagnostics.
* Advance the native state API marker to 3 to detect stale DLLs and workers.
* Retain named NULL refit fields after failed CV, preventing R's partial `$fit`
  matching from incorrectly returning fit_time.
* Correct the optional fork-support probe to parallelly::supportsMulticore;
  the old parallel:: lookup always fell back to independent sessions.
* Compile and test the installed package on R 4.3.3/Linux, including numerical
  comparisons against independent dense updates, cross-backend warm starts,
  rank-1/2/5 checks, and serial/multisession CV. See ALGORITHM_IMPLEMENTATION_0212.md.

# egcar 0.2.11

* Build compact accelerated preparations directly; wide views use a single thin
  SVD rather than first computing an unused full covariance eigendecomposition.
* Allocate dense block-diagonal covariances and full eigenvalue-product caches
  only when dense reference/oracle solvers require them.
* Score wide validation splits through n_validation by rank projections, with
  the same training centering, covariance divisor, ridge and score definition.
* Reuse full-support spectra for loading normalization, including the ridge on
  the complement of a thin basis. Zero-ridge and partial-support cases retain
  the original dense factorization.
* Stream group norms and proximal operations without retaining Wk/Wl edge
  arrays; retain per-edge transformed coefficients only for requested histories.
* Bound the loading-factor cache by both entry count and total estimated bytes
  (64 MiB by default; configurable with `loading_cache_max_bytes`).
* Carry endpoint group states directly along CV paths. Public fitted G/V fields
  and legacy warm starts remain supported.
* Limit EGCAR worker closures to required inputs and omit unused training views
  from EGCAR fold payloads. Shared CV objects still support comparison methods.
* Preserve 0.2.10 future-plan/NSE handling. Reject missing native API markers as
  documented, and repair a legacy test that referenced removed fold fields.
* Add numerical regression tests and a reproducible old/new R benchmark.
  Algebra and extracted C++ kernel checks passed locally; R execution and full
  package benchmarks were unavailable. See OPTIMIZATION_REPORT.md.

# egcar 0.2.10

* Actual fix for the `multicore` crash 0.2.9 failed to resolve (root cause identified from a full, untruncated `traceback()`, not a guess this time): `.egcar_with_workers()` (`R/utils.R`) called `future::plan(if (use_multicore) future::multicore else future::multisession, workers = ...)`. `future::plan()` uses non-standard evaluation on its `strategy` argument (`substitute = TRUE` by default), so it captured the *literal, unevaluated* `if(...)` call rather than the plain function that conditional would resolve to. Its internal `tweak()` then treated the primitive `` `if` `` itself as the strategy being configured and tried to inspect its environment; primitives have none (`environment(`if`)` is `NULL`), producing `Error in ls(envir = env, all.names = TRUE) : invalid 'envir' argument` inside `tweak.function` -- unconditionally, at `plan()` setup, before any CV-fold work is ever reached. This explains why 0.2.9's change (to the `future_lapply()` call inside the fold dispatcher) had no effect: the crash never reached that code at all.
  The conditional is now evaluated into a plain variable (`strategy_fn`) before being passed to `plan()`, so `plan()`'s `substitute()` captures a simple symbol instead of an `if` expression -- exactly how the original, always-working `future::plan(future::multisession, workers = ...)` call was shaped before 0.2.8 introduced the inline conditional. No other `plan()` call sites in the package pass a non-trivial expression as `strategy`.

# egcar 0.2.9

* Attempted fix for a crash introduced by the 0.2.8 `multicore` change (untested against a live R session -- see caveat below): `.egcar_engine()`'s CV-fold dispatcher (`R/globals.R`) calls `future.apply::future_lapply()` on a function (`task`) whose environment was deliberately reassigned to a hand-built "sandbox" environment (`environment(dispatcher) <- e`), used throughout the package so each engine instance gets its own isolated runtime configuration (`BLAS_THREADS`, `EGCAR_PARTIAL_EIGEN`, etc.) without touching `.GlobalEnv`. Left at its default (`future.globals = TRUE`, i.e. auto-detect), `future`/`globals`'s dependency-scanning walked this reassigned environment and crashed under `future::multicore` specifically with `Error in ls(envir = env, all.names = TRUE) : invalid 'envir' argument` (inside `globals`'s internal `tweak.function`) -- even though the identical environment shape had always been tolerated fine under `future::multisession`.
  `future.globals` is now supplied explicitly as a named list (`BLAS_THREADS`, `FUN`, `.egcar_with_threads`) instead, so `future`/`globals` is told exactly what `task` needs rather than having to discover it by walking `e`. This should have no effect on results under either backend: the same values that auto-detection would have found are now just handed over directly.
  Caveat: this could not be tested against a live R/SLURM environment before release. If it does not resolve the crash, `EGCAR_MULTICORE=0` (introduced in 0.2.8) remains available to force `multisession` unconditionally while this is investigated further.

# egcar 0.2.8

* Performance/memory change (not a bug fix): `.egcar_with_workers()` (`R/utils.R`) now uses `future::multicore` instead of `future::multisession` for the parallel CV-fold workers, whenever OS-level forking is actually available (Unix-like systems where `parallel::supportsMulticore()` returns TRUE; Windows and unsupported environments still use `multisession` exactly as before). `multicore` forks the running R process, so each worker shares the parent's already-computed, read-only data (prepared covariance blocks, cached eigenbases, etc.) via copy-on-write, rather than each worker being a fully independent R session with its own complete copy of everything it touches. At large p this is the difference between roughly 1x and roughly `workers`x peak memory for the same 5-fold CV parallelism -- directly relevant to configurations like p_per_block=5000, where 5 independent `multisession` workers each needed to hold their own multi-GB working set.
  Checked before making this change: `future.seed = TRUE` is already set on every `future_lapply()` call in the package, so per-worker random draws stay independent regardless of which backend is used; and every SLURM script in this project already pins `OMP_NUM_THREADS=1`/`MKL_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1` etc., which avoids the classic multithreaded-BLAS-plus-fork deadlock risk that would otherwise make `multicore` riskier than `multisession`.
  Set `EGCAR_MULTICORE=0` to force `multisession` even on a system that supports forking, if ever needed.

# egcar 0.2.7

* Fixed a second, related crash: `egcar_loading_factors()` and `loading_metric_factors()` (`R/helpers.r`) cached their results by `assign()`/`exists()`/`get()` on an environment, keyed by a string built from `paste(selected, collapse = ",")` -- the full comma-joined list of selected row indices. R caps every environment variable name at 10,000 bytes; at p_total in the low thousands (confirmed failing at p_per_block=1000, p_total=3000, where the joined key reached ~13,900 bytes) this crashed with "variable names are limited to 10000 bytes". At p_total=15000 (p_per_block=5000) the same key would have reached ~79,000 bytes, so this would also have blocked the p=5000 configuration even after the 0.2.6 overflow fix, just later in the pipeline (after the ADMM fit, during loading extraction).
  Both functions now key their cache on a new shared helper, `.egcar_cache_key()`, which builds a short, fixed-width surrogate (independent of how many rows are selected) instead of the raw index list. This surrogate is not collision-proof by itself and is not meant to be: the actual `selected` vector and `covariance_ridge` are stored alongside the cached value and checked with `identical()` before a cache hit is trusted, so a collision only costs a harmless recompute, never a silently wrong cached result.

# egcar 0.2.6

* Fixed a 32-bit integer overflow in `egcar_prepare_context()` (`R/helpers.r`): the `project_left`/`lift_left` cost-comparison products (`a * pk * pl`, `pk * pl * b`, etc., where `pk`/`pl` come from `p_list` and `a`/`b` from `lengths()` of the cached spectra) were computed in native R integer arithmetic, which silently overflows to `NA` past roughly p_total = 1,900-2,600 at full rank -- rather than raising an error, R only issues an "NAs produced by integer overflow" warning. The resulting `NA` then reached `if (z$project_left[[e]]) ...` inside `egcar_project()`, failing with "missing value where TRUE/FALSE needed". `pk`, `pl`, `a`, `b` are now coerced to double before these products are formed, which is exact for these magnitudes and eliminates the overflow. Confirmed via the affected p_per_block=5000 (p_total=15000) configuration.
* No other integer-overflow-prone triple products of dimension-like quantities were found elsewhere in the R sources; the compiled backend (`src/egcar_native.cpp`) uses Armadillo's 64-bit `arma::uword` throughout and is not subject to the same risk.

# egcar 0.2.5

* Harden experiment configuration compatibility: retention-only fields missing from older config objects are filled with the current memory-first defaults before validation.
* Namespace-qualify the installed small all-methods example so stale functions in `.GlobalEnv` cannot mask `egcar::egcar_experiment_config()` or `egcar::run_egcar_experiments()`.
* The small diagnostic explicitly retains fold-level CV diagnostics.

# egcar 0.2.4

* Hardened the optimized/native solver boundary against legacy nested Rcpp conversions that may return an edge matrix as a dimensionless numeric vector. Solver state is now normalized to the known edge dimensions before any row/column reductions.
* `egcar_group_norms()` and `egcar_view_copies()` now validate/restore edge-matrix dimensions defensively. This is a compatibility safeguard only; numerical values, column-major ordering, ADMM equations, penalties, CV rules, and estimators are unchanged.
* The small all-methods and matrix-interface examples now require 0.2.4.

# egcar 0.2.3

* Memory-first experiment defaults: full fitted objects, fold-level CV tables, full loading-plot data, and per-configuration loading PDFs are no longer retained/generated unless explicitly requested.
* Save compact per-configuration loading archives with xz compression so loading heatmaps can be reconstructed without complete solver objects.
* Drop unused fold indices/training means and dispatch fold objects directly to workers to reduce future-global export pressure.
* Compact final experiment ADMM results: rate fits and selected refits keep C_hat and diagnostics but discard warm-start auxiliary state.
* For group ADMM, compute active rows directly from endpoint copies and avoid constructing view-wide G/V matrices unless state retention is requested.
* Compact SGCA caches by retaining rank-r U plus scalar initializer diagnostics rather than Pi/H/Gamma; remove an unused cached U0 matrix.
* Suppress complete external benchmark fit objects in the experiment by default; public comparison calls retain their previous behavior.
* Reduce the default loading-factor cache from 128 to 4 entries, default `keep_full_C=FALSE` because edge blocks already represent the operator, construct projector differences on demand, stop writing duplicate plot-data CSVs, and use xz-compressed checkpoints.
* Numerical estimators, tuning grids, folds, seeds, validation losses and parameter-selection rules are unchanged.

# egcar 0.2.2

* Explicitly construct each native state block as an R numeric matrix, including
  rectangular, one-row, one-column, and 1-by-1 blocks.
* Check the native matrix API version and returned state shapes before computing
  row norms or extracting loadings; do not silently reshape dimensionless arrays.
* Apply the same serialization fix to the bundled standalone runner and use a
  new native-cache filename/option there.
* Add installed matrix-interface, rank 1/2/5, warm-start, rate and CV regression
  checks, plus focused testthat coverage.
* Leave numerical solver updates, penalties, grids, seeds, benchmarks, oracles,
  experiment defaults and parallel-worker allocation unchanged.
* R/Rcpp execution of these new checks has not been performed in the patching
  environment. See inst/doc/MATRIX_INTERFACE_FIX_022.md for validation scope.

# egcar 0.2.1

* Speed-only revision of the compiled EGCAR L11/L21 backend.
* Adapted mapped Eigen products for small matrices after reviewing EfficientCCA.
* Kept Armadillo kernels for tiny and large products; disabled Eigen threading.
* Removed copies of read-only double-valued native context matrices.
* Reused projection, lifting and group-norm buffers; fused proximal/dual/residual loops.
* Added RcppEigen as a build dependency, without requiring EfficientCCA or SMUT.
* Preserved every R implementation, original controls, other methods, CV and oracles.
* Added a one-configuration, ten-method check using five shared folds/workers.
* Added native old/new comparisons, standalone timing checks and R-level tests.

# egcar 0.2.0

* Separate explicit EGCAR_L11 and EGCAR_L21 fit/CV/rate entry points.
* Strictly positive EGCAR CV grids, including the full and smoke experiments.
* Removed the redundant tied experiment and legacy result names.
* Bundled the four required SGCA initializer functions and retained MIT notices.
* Retained the local experiment's SGCA/RGCCA/SGCCA/MultiCCA CV routines.
* Added a shared-config full experiment runner and matching standalone source.
* Added explicit CV-failure diagnostics, positive-grid and initializer tests.
* Fixed the stats::toeplitz namespace import missing in the prior package.
* Retained the accelerated C++ core, dense references, and original oracles.

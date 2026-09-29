\name{egcar_control}
\alias{egcar_control}
\alias{benchmark_control}
\title{Controls for EGCAR and comparison solvers}
\description{Validated per-call numerical controls. They do not alter package defaults, global variables, or future plans.}
\usage{
egcar_control(
    backend = c("cpp", "R", "reference"), max_iter = 2000L,
    max_iter_cv = 1000L, abs_tol = 1e-5, rel_tol = 1e-4,
    mu = 1, adaptive_mu = TRUE, balance_ratio = 10, scale_factor = 2,
    adapt_every = 10L, check_every = 1L, check_every_cv = 5L,
    row_threshold = 1e-4, covariance_ridge = 1e-4,
    group_zero_tol = 1e-8, entry_zero_tol = 1e-10,
    partial_eigen = TRUE, partial_eigen_min = 2L,
    loading_cache_max = 4L, keep_history = FALSE,
    keep_full_C = FALSE, blas_threads = 1L, verbose = FALSE,
  loading_cache_max_bytes = 64 * 1024^2, compact_state = TRUE,
  cv_require_convergence = TRUE)

benchmark_control(
    sgca_eta = 0.001, sgca_ridge = 0, sgca_init_tol = 5e-3,
    sgca_init_max_iter = 1000L, sgca_tgd_tol = 1e-6,
    sgca_tgd_max_iter = 15000L, fast_sgca_initializer = TRUE,
    rgcca_scheme = "factorial", rgcca_tol = 1e-8, rgcca_max_iter = 1000L,
    multicca_niter = 25L, multicca_backend = c("gram", "PMA"),
    align_signs = TRUE, sgca_time_limit = 6 * 60 * 60,
    sgca_stopping = c("fixed_iterations", "absolute_change"),
    sgca_cv_score = c("paper", "common_loss"),
    sgca_cv_require_convergence = FALSE)
}
\arguments{
\item{backend}{Compiled \code{"cpp"}, accelerated \code{"R"}, or dense \code{"reference"}. All use the same statistical objectives. A source installation still needs a compiler even when selecting the R backend.}
\item{max_iter, max_iter_cv}{ADMM iteration caps for final fits and training folds.}
\item{abs_tol, rel_tol}{Absolute and relative full-space ADMM residual tolerances.}
\item{mu}{Initial augmentation parameter; the group-only edge shift is twice this value.}
\item{adaptive_mu}{Whether to use the supplied residual-balancing heuristic.}
\item{balance_ratio, scale_factor, adapt_every}{Residual balancing ratio, multiplicative change, and adaptation frequency. Adaptation occurs only on iterations that also check residuals.}
\item{check_every, check_every_cv}{Residual-check frequencies for final fits and CV folds. The first and last iterations are always checked.}
\item{row_threshold}{Row-norm threshold used for localized loading extraction.}
\item{covariance_ridge}{Ridge multiplier times the mean selected marginal variance. It regularizes loading normalization, not the regression loss.}
\item{group_zero_tol, entry_zero_tol}{Numerical output cutoffs inherited from the supplied implementation; not additional regularization penalties.}
\item{partial_eigen, partial_eigen_min}{Enable matrix-free RSpectra loading eigensolves and set their minimum selected dimension. The default computes only the leading rank eigenpairs when rank is smaller than the selected dimension. Failure never triggers an automatic dense fallback. FALSE explicitly enables dense diagnostic extraction.}
\item{loading_cache_max}{Maximum number of cached support-specific loading factorizations.}
\item{loading_cache_max_bytes}{Maximum total estimated bytes in the loading-factor cache. Defaults to 64 MiB; zero disables caching. This does not limit active-factorization or solver memory.}
\item{cv_require_convergence}{Require convergence and finite scores on every EGCAR CV fold. SGCA uses its separate sgca_cv_require_convergence control.}
\item{compact_state}{Retain compressed Hk/Hl/a group warm-start state. FALSE reconstructs legacy view-wide G/V output; optimized iterations still use compressed state.}
\item{keep_history}{Record the group-only ADMM objective history. The entrywise solvers do not record histories.}
\item{keep_full_C}{Retain the assembled operator with valid loading results. Edge blocks are always retained.}
\item{blas_threads}{Optional scoped BLAS/OpenMP thread limit when RhpcBLASctl is installed. NULL leaves thread settings untouched.}
\item{verbose}{Print solver progress.}
\item{sgca_eta, sgca_ridge}{SGCA paper step parameter (the gradient-half multiplier is 2*eta) and optional metric ridge. Ridge defaults to zero.}
\item{sgca_init_tol, sgca_init_max_iter}{Absolute Frobenius initializer change tolerance and outer iteration cap; defaults 0.005 and 1000 from the authors' R controls. Inf removes the cap.}
\item{sgca_tgd_tol, sgca_tgd_max_iter}{Absolute Frobenius TGD iterate-change threshold and number of updates. Defaults 1e-6 and 15000. The threshold is diagnostic only in fixed_iterations mode. Inf is allowed only with absolute_change stopping.}
\item{sgca_stopping}{Fixed paper iterations, or optional absolute-change early stopping.}
\item{sgca_cv_score}{Paper covariance-trace score, or the common EGCAR validation loss.}
\item{sgca_cv_require_convergence}{Require both numerical change tests on all SGCA folds and refit. FALSE permits finite capped initialization and completed fixed-step TGD, with diagnostics preserved.}
\item{sgca_time_limit}{Elapsed seconds for one complete SGCA CV grid and final refit, including initialization. Defaults to six hours. A supervised R process can terminate compiled operations. Shared data preparation precedes the budget. Inf disables supervision. Startup and process cleanup can add overhead.}
\item{fast_sgca_initializer}{Use cached matrix algebra for the bundled SGCA initializer; FALSE selects its direct reference implementation.}
\item{rgcca_scheme, rgcca_tol, rgcca_max_iter}{RGCCA-family package solver controls. The supplied benchmark uses factorial scheme.}
\item{multicca_niter, multicca_backend}{MultiCCA iteration limit and cached-Gram or public-PMA implementation.}
\item{align_signs}{Resolve external block signs using training data only before computing the common signed score.}
}
\value{A named control list. Benchmark controls affect only comparison methods.}
\examples{
egcar_control()
benchmark_control()
}

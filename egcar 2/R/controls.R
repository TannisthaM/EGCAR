# Validated per-call EGCAR and comparison-method controls.

#' Controls for EGCAR and comparison solvers
#'
#' @description Validated per-call numerical controls. They do not alter package defaults, global variables, or future plans.
#' @param backend Compiled \code{"cpp"}, accelerated \code{"R"}, or dense \code{"reference"}. All use the same statistical objectives. A source installation still needs a compiler even when selecting the R backend.
#' @param max_iter,max_iter_cv ADMM iteration caps for final fits and training folds.
#' @param abs_tol,rel_tol Absolute and relative full-space ADMM residual tolerances.
#' @param mu Initial augmentation parameter; the group-only edge shift is twice this value.
#' @param adaptive_mu Whether to use the supplied residual-balancing heuristic.
#' @param balance_ratio,scale_factor,adapt_every Residual balancing ratio, multiplicative change, and adaptation frequency. Adaptation occurs only on iterations that also check residuals.
#' @param check_every,check_every_cv Residual-check frequencies for final fits and CV folds. The first and last iterations are always checked.
#' @param row_threshold Row-norm threshold used for localized loading extraction.
#' @param covariance_ridge Ridge multiplier times the mean selected marginal variance. It regularizes loading normalization, not the regression loss.
#' @param group_zero_tol,entry_zero_tol Numerical output cutoffs inherited from the supplied implementation; not additional regularization penalties.
#' @param partial_eigen,partial_eigen_min Enable checked RSpectra loading eigensolves and set their minimum selected dimension. Dense fallback is used when unavailable or a check fails.
#' @param loading_cache_max Maximum number of cached support-specific loading factorizations.
#' @param keep_history Record the group-only ADMM objective history. The entrywise solvers do not record histories.
#' @param keep_full_C Retain the assembled operator with valid loading results. Edge blocks are always retained.
#' @param blas_threads Optional scoped BLAS/OpenMP thread limit when RhpcBLASctl is installed. NULL leaves thread settings untouched.
#' @param verbose Print solver progress.
#' @param sgca_eta,sgca_ridge SGCA gradient step parameter and initializer metric ridge.
#' @param sgca_init_tol,sgca_init_max_iter SGCA initializer tolerance and iteration limit.
#' @param sgca_tgd_tol,sgca_tgd_max_iter Penalized SGCA gradient-mapping tolerance and iteration limit.
#' @param fast_sgca_initializer Use cached matrix algebra for the bundled SGCA initializer; FALSE selects its direct reference implementation.
#' @param rgcca_scheme,rgcca_tol,rgcca_max_iter RGCCA-family package solver controls. The supplied benchmark uses factorial scheme.
#' @param multicca_niter,multicca_backend MultiCCA iteration limit and cached-Gram or public-PMA implementation.
#' @param align_signs Resolve external block signs using training data only before computing the common signed score.
#' @return A named control list. Benchmark controls affect only comparison methods.
#' @rdname egcar_control
#' @examples
#' egcar_control()
#' benchmark_control()
#' @export
egcar_control <- function(
    backend = c("cpp", "R", "reference"), max_iter = 2000L,
    max_iter_cv = 1000L, abs_tol = 1e-5, rel_tol = 1e-4,
    mu = 1, adaptive_mu = TRUE, balance_ratio = 10, scale_factor = 2,
    adapt_every = 10L, check_every = 1L, check_every_cv = 5L,
    row_threshold = 1e-4, covariance_ridge = 1e-4,
    group_zero_tol = 1e-8, entry_zero_tol = 1e-10,
    partial_eigen = TRUE, partial_eigen_min = 128L,
    loading_cache_max = 128L, keep_history = FALSE,
    keep_full_C = TRUE, blas_threads = 1L, verbose = FALSE) {
  if (length(backend) == 1L && identical(backend, "r")) backend <- "R"
  backend <- match.arg(backend)
  for (nm in c("max_iter", "max_iter_cv", "adapt_every", "check_every",
               "check_every_cv", "partial_eigen_min"))
    .egcar_scalar(get(nm), nm, 1, integer = TRUE)
  .egcar_scalar(loading_cache_max, "loading_cache_max", 0, integer = TRUE)
  for (nm in c("abs_tol", "rel_tol", "row_threshold", "covariance_ridge",
               "group_zero_tol", "entry_zero_tol"))
    .egcar_scalar(get(nm), nm)
  for (nm in c("mu", "balance_ratio")) .egcar_scalar(get(nm), nm, strict = TRUE)
  .egcar_scalar(scale_factor, "scale_factor", 1, strict = TRUE)
  for (nm in c("adaptive_mu", "partial_eigen", "keep_history", "keep_full_C", "verbose"))
    .egcar_flag(get(nm), nm)
  if (!is.null(blas_threads)) .egcar_scalar(blas_threads, "blas_threads", 1, integer = TRUE)
  if (backend == "reference" && keep_history)
    message("The reference l11 solver does not record objective histories.")
  values <- mget(names(formals(egcar_control)), envir = environment())
  for (nm in c("max_iter", "max_iter_cv", "adapt_every", "check_every",
               "check_every_cv", "partial_eigen_min", "loading_cache_max"))
    values[[nm]] <- as.integer(values[[nm]])
  class(values) <- "egcar_control"
  values
}

#' @rdname egcar_control
#' @export
benchmark_control <- function(
    sgca_eta = 0.001, sgca_ridge = 1e-6, sgca_init_tol = 5e-3,
    sgca_init_max_iter = 1000L, sgca_tgd_tol = 1e-6,
    sgca_tgd_max_iter = 15000L, fast_sgca_initializer = TRUE,
    rgcca_scheme = "factorial", rgcca_tol = 1e-8, rgcca_max_iter = 1000L,
    multicca_niter = 25L, multicca_backend = c("gram", "PMA"),
    align_signs = TRUE) {
  multicca_backend <- match.arg(multicca_backend)
  rgcca_scheme <- match.arg(rgcca_scheme, c("factorial", "centroid", "horst"))
  for (nm in c("sgca_eta", "sgca_init_tol", "sgca_tgd_tol", "rgcca_tol"))
    .egcar_scalar(get(nm), nm, strict = TRUE)
  .egcar_scalar(sgca_ridge, "sgca_ridge")
  for (nm in c("sgca_init_max_iter", "sgca_tgd_max_iter", "rgcca_max_iter", "multicca_niter"))
    .egcar_scalar(get(nm), nm, 1, integer = TRUE)
  .egcar_flag(fast_sgca_initializer, "fast_sgca_initializer")
  .egcar_flag(align_signs, "align_signs")
  values <- mget(names(formals(benchmark_control)), envir = environment())
  class(values) <- "egcar_benchmark_control"
  values
}

.egcar_as_control <- function(x) {
  if (!is.list(x)) stop("control must be an egcar_control() object or a named list.")
  if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))
    stop("control must have unique, nonempty names.")
  unknown <- setdiff(names(x), names(formals(egcar_control)))
  if (length(unknown)) stop("Unknown control: ", paste(unknown, collapse = ", "))
  do.call(egcar_control, x)
}

.egcar_as_benchmark_control <- function(x) {
  if (!is.list(x)) stop("benchmarks must be a benchmark_control() object or named list.")
  if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))
    stop("benchmarks must have unique, nonempty names.")
  unknown <- setdiff(names(x), names(formals(benchmark_control)))
  if (length(unknown)) stop("Unknown benchmark control: ", paste(unknown, collapse = ", "))
  do.call(benchmark_control, x)
}

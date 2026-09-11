# Shared defaults for standalone and package experiment runners.
egcar_positive_cv_grid <- function(x, name = "lambda") {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop(name, " must contain finite, strictly positive coefficients; zero is not a CV candidate.",
         call. = FALSE)
  }
  unique(as.numeric(x))
}

validate_egcar_experiment_config <- function(config) {
  if (!is.list(config) || is.null(names(config)) || anyDuplicated(names(config)))
    stop("config must be a uniquely named configuration list.")
  for (nm in c("p_list", "n_grid", "rank_grid")) {
    x <- config[[nm]]
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x != floor(x) | x < 1))
      stop(nm, " must contain positive finite integers.")
    config[[nm]] <- as.integer(x)
  }
  if (length(config$p_list) < 2L) stop("At least two views are required.")
  if (anyDuplicated(config$n_grid) || anyDuplicated(config$rank_grid))
    stop("n_grid and rank_grid must not contain duplicates.")
  config$rho_e_cv_grid <- egcar_positive_cv_grid(config$rho_e_cv_grid, "rho_e_cv_grid")
  config$lambda_g_cv_grid <- egcar_positive_cv_grid(config$lambda_g_cv_grid, "lambda_g_cv_grid")
  positive_integer <- c("active_per_view", "n_folds", "max_iter_cv", "max_iter_final",
    "check_every_admm", "sgca_max_iter_init", "sgca_max_iter_tgd", "rgcca_max_iter",
    "multicca_niter", "oracle1_max_iter", "loading_methods_per_page")
  for (nm in positive_integer) {
    v <- config[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 || v != floor(v))
      stop(nm, " must be a positive integer.")
    config[[nm]] <- as.integer(v)
  }
  for (nm in c("master_seed", "loading_factor_cache_max")) {
    v <- config[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 || v != floor(v))
      stop(nm, " must be a nonnegative integer.")
  }
  if (config$n_folds < 2L || any(config$n_grid <= config$n_folds))
    stop("Use at least two folds and n greater than n_folds.")
  if (any(config$rank_grid > pmin(config$active_per_view, min(config$p_list))))
    stop("Each rank must be <= active_per_view and every view dimension.")
  if (length(config$toeplitz_rho) != 1L && length(config$toeplitz_rho) != length(config$p_list))
    stop("toeplitz_rho must have length one or one entry per view.")
  if (any(!is.finite(config$toeplitz_rho)) || any(abs(config$toeplitz_rho) >= 1))
    stop("Toeplitz correlations must lie strictly between -1 and 1.")
  if (length(config$signal) != 1L || !is.finite(config$signal) || config$signal <= 0 || config$signal >= 1)
    stop("signal must lie strictly between zero and one.")
  for (nm in c("mu_z", "mu_g")) {
    x <- config[[nm]]
    if (length(x) != 1L || !is.finite(x) || x <= 0) stop(nm, " must be positive.")
  }
  for (nm in c("abs_tol", "rel_tol", "row_threshold", "covariance_ridge",
               "group_zero_tol", "entry_zero_tol", "rate_c_e", "rate_c_g")) {
    x <- config[[nm]]
    if (length(x) != 1L || !is.finite(x) || x < 0) stop(nm, " must be nonnegative.")
  }
  for (nm in c("adaptive_mu", "run_external_benchmarks", "stop_if_benchmark_packages_missing",
               "align_external_block_signs", "save_fits", "save_cv_fold_results",
               "save_loading_data", "save_compact_loadings", "retain_benchmark_fits",
               "make_plots", "fast_sgca_initializer", "make_loading_plots")) {
    if (!is.logical(config[[nm]]) || length(config[[nm]]) != 1L || is.na(config[[nm]]))
      stop(nm, " must be TRUE or FALSE.")
  }
  if (!config$multicca_backend %in% c("gram", "PMA")) stop("Invalid multicca_backend.")
  class(config) <- "egcar_experiment_config"
  config
}

egcar_experiment_config <- function(
    p_list = c(15L, 15L, 15L),
    n_grid = c(30L, 45L, 100L, 1000L, 5000L, 10000L),
    rank_grid = c(1L, 2L, 5L), ...) {
  P_LIST <- p_list
  N_GRID <- n_grid
  RANK_GRID <- rank_grid
  ACTIVE_PER_VIEW <- 5L
  TOEPLITZ_RHO <- c(0.5, 0.7, 0.9)
  SIGNAL <- 0.8
  MASTER_SEED <- 20260907L

  # Direct penalty grids used in cross-validation.  These are the actual values
  # of rho_e and lambda_g; they are not multiplied by n-dependent rates.
  RHO_E_CV_GRID <- 10^seq(-5, 4)
  LAMBDA_G_CV_GRID <- 10^seq(-5, 4)

  # Multipliers used only by the separate rate-scaled benchmark curves below.
  # They do not enter cross-validation.
  RATE_C_E <- 1
  RATE_C_G <- 1
  N_FOLDS <- 5L

  MU_Z <- 1
  MU_G <- 1
  MAX_ITER_CV <- 1000L
  MAX_ITER_FINAL <- 2000L
  ABS_TOL <- 1e-5
  REL_TOL <- 1e-4
  ADAPTIVE_MU <- TRUE
  ROW_THRESHOLD <- 1e-4
  COVARIANCE_RIDGE <- 1e-4
  GROUP_ZERO_TOL <- 1e-8
  ENTRY_ZERO_TOL <- 1e-10

  # ADMM residuals/convergence/adaptive-mu are only recomputed every
  # CHECK_EVERY_ADMM iterations (plus the first and last), instead of every
  # single iteration, to cut R-level overhead across the very large number of
  # ADMM solves performed during cross-validation.
  CHECK_EVERY_ADMM <- 5L

  # Legacy coarse-to-fine settings are retained for interface compatibility.
  # The replacement CV has only a lambda_g dimension and evaluates its full
  # grid directly, so the two-dimensional coarse/refinement branch is unused.
  COARSE_STEP_CV <- 2L
  REFINE_WINDOW_CV <- 1L

  # External benchmark methods are enabled by default.  Set this to FALSE only
  # for a quick diagnostic run of the new estimator.
  RUN_EXTERNAL_BENCHMARKS <- TRUE
  STOP_IF_BENCHMARK_PACKAGES_MISSING <- TRUE
  # SGCA searches the full Cartesian product of these three grids. k counts
  # nonzero rows across ALL views, not per view. Invalid k < rank or k > p are
  # dropped separately for each rank; no population support is used to select k.
  SGCA_K_GRID <- sort(unique(c(5L, 10L, 15L, 20L, 30L, sum(P_LIST))))
  SGCA_RHO_GRID <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1)
  SGCA_LAMBDA_GRID <- 10^seq(-5, 4)
  SGCA_ETA <- 0.001
  SGCA_RIDGE_B <- 1e-6
  SGCA_INIT_TOL <- 5e-3
  SGCA_MAX_ITER_INIT <- 1000L
  SGCA_TGD_TOL <- 1e-6
  SGCA_MAX_ITER_TGD <- 15000L

  RGCCA_TAU_GRID <- c(1e-6, 1e-3, 0.1, 0.25, 0.5, 0.75, 1)
  RGCCA_SCHEME <- "factorial"
  RGCCA_TOL <- 1e-8
  RGCCA_MAX_ITER <- 1000L
  # These are common scalar bounds across views, with valid package ranges.
  # Larger SGCCA sparsity / MultiCCA L1 bounds mean LESS regularization.
  SGCCA_SPARSITY_GRID <- seq(max(1 / sqrt(P_LIST)), 1, length.out = 10L)
  MULTICCA_L1_GRID <- seq(1, min(sqrt(P_LIST)), length.out = 10L)
  MULTICCA_NITER <- 25L
  ALIGN_EXTERNAL_BLOCK_SIGNS <- TRUE

  ORACLE1_MAX_ITER <- MAX_ITER_FINAL
  # Memory/storage defaults: retain scientific summaries and compact loadings,
  # not complete solver/package objects or fold-level diagnostics.
  SAVE_FITS <- FALSE
  SAVE_CV_FOLD_RESULTS <- FALSE
  SAVE_LOADING_DATA <- FALSE
  SAVE_COMPACT_LOADINGS <- TRUE
  RETAIN_BENCHMARK_FITS <- FALSE
  MAKE_PLOTS <- TRUE

  # Matrix-only accelerations. Reference backends remain available for checks.
  FAST_SGCA_INITIALIZER <- TRUE
  MULTICCA_BACKEND <- "gram"  # "gram" (fast) or "PMA" (original package call)
  LOADING_FACTOR_CACHE_MAX <- 4L  # small bounded support-specific whitening cache

  # Separate loading PDFs for EVERY completed (rep, rank, n), by default.
  # NULL selects all; e.g. LOADING_PLOT_N <- c(100L, 10000L) limits PDF output.
  # Per-configuration loading PDFs are large and can be regenerated from
  # compact_loadings/, so they are off in the memory-first default.
  MAKE_LOADING_PLOTS <- FALSE
  LOADING_PLOT_N <- NULL
  LOADING_PLOT_RANKS <- NULL
  LOADING_PLOT_REPS <- NULL
  LOADING_METHODS_PER_PAGE <- 3L  # truth is repeated on each comparison page

  keys <- c("P_LIST", "N_GRID", "RANK_GRID", "ACTIVE_PER_VIEW", "TOEPLITZ_RHO", "SIGNAL", "MASTER_SEED", "RHO_E_CV_GRID", "LAMBDA_G_CV_GRID", "RATE_C_E", "RATE_C_G", "N_FOLDS", "MU_Z", "MU_G", "MAX_ITER_CV", "MAX_ITER_FINAL", "ABS_TOL", "REL_TOL", "ADAPTIVE_MU", "ROW_THRESHOLD", "COVARIANCE_RIDGE", "GROUP_ZERO_TOL", "ENTRY_ZERO_TOL", "CHECK_EVERY_ADMM", "COARSE_STEP_CV", "REFINE_WINDOW_CV", "RUN_EXTERNAL_BENCHMARKS", "STOP_IF_BENCHMARK_PACKAGES_MISSING", "SGCA_K_GRID", "SGCA_RHO_GRID", "SGCA_LAMBDA_GRID", "SGCA_ETA", "SGCA_RIDGE_B", "SGCA_INIT_TOL", "SGCA_MAX_ITER_INIT", "SGCA_TGD_TOL", "SGCA_MAX_ITER_TGD", "RGCCA_TAU_GRID", "RGCCA_SCHEME", "RGCCA_TOL", "RGCCA_MAX_ITER", "SGCCA_SPARSITY_GRID", "MULTICCA_L1_GRID", "MULTICCA_NITER", "ALIGN_EXTERNAL_BLOCK_SIGNS", "ORACLE1_MAX_ITER", "SAVE_FITS", "SAVE_CV_FOLD_RESULTS", "SAVE_LOADING_DATA", "SAVE_COMPACT_LOADINGS", "RETAIN_BENCHMARK_FITS", "MAKE_PLOTS", "FAST_SGCA_INITIALIZER", "MULTICCA_BACKEND", "LOADING_FACTOR_CACHE_MAX", "MAKE_LOADING_PLOTS", "LOADING_PLOT_N", "LOADING_PLOT_RANKS", "LOADING_PLOT_REPS", "LOADING_METHODS_PER_PAGE")
  out <- mget(keys, envir = environment(), inherits = FALSE)
  names(out) <- tolower(names(out))
  changes <- list(...)
  if (length(changes)) {
    if (is.null(names(changes)) || any(!nzchar(names(changes))) || anyDuplicated(names(changes)))
      stop("Configuration overrides must have unique, nonempty names.")
    unknown <- setdiff(names(changes), names(out))
    if (length(unknown)) stop("Unknown configuration field: ", paste(unknown, collapse = ", "))
    out[names(changes)] <- changes
  }
  validate_egcar_experiment_config(out)
}

egcar_smoke_config <- function(config) {
  config$n_grid <- 45L
  config$rank_grid <- c(1L, 2L, 5L)
  config$rho_e_cv_grid <- c(0.01, 0.1)
  config$lambda_g_cv_grid <- c(0.01, 0.1)
  config$sgca_k_grid <- c(10L, 20L)
  config$sgca_rho_grid <- c(0.01, 0.1)
  config$sgca_lambda_grid <- c(0.01, 1)
  config$rgcca_tau_grid <- c(0.1, 1)
  config$sgcca_sparsity_grid <- c(0.5, 1)
  config$multicca_l1_grid <- c(2, min(sqrt(config$p_list)))
  config$max_iter_cv <- 30L; config$max_iter_final <- 60L
  config$oracle1_max_iter <- 60L; config$sgca_max_iter_init <- 30L
  config$sgca_max_iter_tgd <- 60L; config$rgcca_max_iter <- 60L
  config$multicca_niter <- 10L
  message("SMOKE TEST: short iteration limits and tiny grids; not performance results.")
  validate_egcar_experiment_config(config)
}

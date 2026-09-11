#!/usr/bin/env Rscript
# Standalone accelerated EGCAR experiment, revision 0.2.0.
# Exactly four EGCAR curves: L11-CV, L11-rate, L21-CV, L21-rate.
# Positive CV grids only; no combined/tied EGCAR benchmark.
# SGCA initialization is bundled with its MIT notice; its local CV/TGD are retained.
# Other comparisons: Oracle1-population, Oracle2-support, SGCA, RGCCA, SGCCA, MultiCCA.
# Required for all comparisons/graphics: RGCCA, PMA, ggplot2.
# Fast native backend: Rcpp, RcppArmadillo, and a C++14 compiler.
# Optional: future, future.apply, RSpectra, RhpcBLASctl.
# Rscript run_EGCAR_local.R OUTPUT_DIR WORKERS N_REPS
# RStudio: source("run_EGCAR_local.R"); run_local_egcar_experiments(...)
# Full defaults reproduce the original simulation design, not a reduced demo.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

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
  retention_defaults <- list(
    save_fits = FALSE, save_cv_fold_results = FALSE,
    save_loading_data = FALSE, save_compact_loadings = TRUE,
    retain_benchmark_fits = FALSE
  )
  for (nm in names(retention_defaults)) {
    if (is.null(config[[nm]])) config[[nm]] <- retention_defaults[[nm]]
  }
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

run_local_egcar_experiments <- function(
    output_dir = "egcar_local_outputs", workers = 1L, n_reps = 1L,
    config = egcar_experiment_config(),
    backend = Sys.getenv("EGCAR_BACKEND", "auto"),
    smoke_test = identical(Sys.getenv("EGCAR_SMOKE_TEST", "0"), "1"),
    definitions_only = FALSE) {
  config <- validate_egcar_experiment_config(config)
  for (nm in c("workers", "n_reps")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x))
      stop(nm, " must be a positive integer.")
  }
  backend <- match.arg(tolower(backend), c("auto", "cpp", "r", "reference"))
  if (smoke_test) { config <- egcar_smoke_config(config); n_reps <- 1L }
  list2env(setNames(unclass(config), toupper(names(config))), envir = environment())
  N_REP <- as.integer(n_reps); SMOKE_TEST <- isTRUE(smoke_test)
  OUT_DIR <- output_dir
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUT_DIR)) stop("Could not create output_dir.")
  OUT_DIR <- normalizePath(OUT_DIR, mustWork = TRUE)
  old_options <- options(stringsAsFactors = FALSE)
  on.exit(options(old_options), add = TRUE)
  N_WORKERS_REQUESTED <- as.integer(workers)
  CV_WORKERS <- as.integer(min(workers, N_FOLDS))
  PARALLEL_CV <- CV_WORKERS > 1L && requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)
  if (PARALLEL_CV) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    old_max <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_max), add = TRUE)
    options(future.globals.maxSize = max(4 * 1024^3, old_max %||% 0))
    future::plan(future::multisession, workers = CV_WORKERS)
  } else {
    if (CV_WORKERS > 1L) message("Parallel packages missing: all CV methods use one worker.")
    CV_WORKERS <- 1L
  }
  N_WORKERS <- CV_WORKERS
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    old_omp <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
    on.exit({
      if (!is.null(old_blas)) try(RhpcBLASctl::blas_set_num_threads(old_blas), silent = TRUE)
      if (!is.null(old_omp)) try(RhpcBLASctl::omp_set_num_threads(old_omp), silent = TRUE)
    }, add = TRUE)
    RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
  }
  EXTERNAL_BENCHMARK_METHODS <- c("SGCA", "RGCCA", "SGCCA", "MultiCCA")
  METHOD_ORDER <- c("Oracle1-population", "Oracle2-support", "MultiCCA",
    "EGCAR-L11-rate", "EGCAR-L11-CV", "EGCAR-L21-rate", "EGCAR-L21-CV",
    "SGCA", "RGCCA", "SGCCA")
  METHOD_COLORS <- setNames(grDevices::hcl.colors(length(METHOD_ORDER), palette = "Dark 3"), METHOD_ORDER)
  CV_METHODS <- c("EGCAR-L11-CV", "EGCAR-L21-CV", EXTERNAL_BENCHMARK_METHODS)
  CV_WORKERS_BY_METHOD <- setNames(rep.int(CV_WORKERS, length(CV_METHODS)), CV_METHODS)
  CV_WORKER_ALLOCATION <- data.frame(method = CV_METHODS,
    requested_workers = rep.int(N_WORKERS_REQUESTED, length(CV_METHODS)),
    allocated_cv_workers = unname(CV_WORKERS_BY_METHOD),
    folds = rep.int(N_FOLDS, length(CV_METHODS)),
    parallel_cv = rep.int(PARALLEL_CV, length(CV_METHODS)),
    enabled = !CV_METHODS %in% EXTERNAL_BENCHMARK_METHODS | RUN_EXTERNAL_BENCHMARKS,
    stringsAsFactors = FALSE)
  utils::write.csv(CV_WORKER_ALLOCATION, file.path(OUT_DIR, "cv_worker_allocation.csv"), row.names = FALSE)
  if (RUN_EXTERNAL_BENCHMARKS) {
    required <- c("RGCCA", "PMA")
    missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing)) {
      msg <- paste("Missing comparison dependencies:", paste(missing, collapse = ", "))
      if (STOP_IF_BENCHMARK_PACKAGES_MISSING) stop(msg) else warning(msg)
    }
  }
  cat(sprintf("Output: %s\nCV workers per method: %d of %d requested; folds=%d; repetitions=%d\n",
    OUT_DIR, CV_WORKERS, N_WORKERS_REQUESTED, N_FOLDS, N_REP))

  set_blas_threads_one <- function() {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
    }
    invisible(NULL)
  }

  # =============================================================================
  # 1. Matrix and block utilities
  # =============================================================================

  symmetrize <- function(A) {
    (A + t(A)) / 2
  }

  frob <- function(A) {
    sqrt(sum(A * A))
  }

  row_l2 <- function(A) {
    sqrt(rowSums(A * A))
  }

  soft_threshold <- function(A, threshold) {
    if (threshold <= 0) return(A)
    sign(A) * pmax(abs(A) - threshold, 0)
  }

  row_group_threshold <- function(A, threshold) {
    if (threshold <= 0) return(A)
    nr <- row_l2(A)
    mult <- pmax(0, 1 - threshold / pmax(nr, .Machine$double.eps))
    A * mult
  }

  matrix_power_psd <- function(A, power, ridge = 0, eig_floor = 1e-10) {
    ev <- eigen(symmetrize(A), symmetric = TRUE)
    values <- pmax(ev$values + ridge, eig_floor)
    tcrossprod(sweep(ev$vectors, 2L, values^power, "*"), ev$vectors)
  }

  block_diag <- function(blocks) {
    nr <- sum(vapply(blocks, nrow, integer(1L)))
    nc <- sum(vapply(blocks, ncol, integer(1L)))
    out <- matrix(0, nr, nc)
    r0 <- 1L
    c0 <- 1L
    for (B in blocks) {
      rr <- r0:(r0 + nrow(B) - 1L)
      cc <- c0:(c0 + ncol(B) - 1L)
      out[rr, cc] <- B
      r0 <- max(rr) + 1L
      c0 <- max(cc) + 1L
    }
    out
  }

  make_block_indices <- function(p_list) {
    starts <- cumsum(c(1L, head(p_list, -1L)))
    lapply(seq_along(p_list), function(k) {
      seq.int(starts[[k]], length.out = p_list[[k]])
    })
  }

  edge_key <- function(k, l) {
    if (k > l) {
      tmp <- k
      k <- l
      l <- tmp
    }
    paste0(k, "_", l)
  }

  make_edge_table <- function(p_list) {
    K <- length(p_list)
    ij <- which(upper.tri(matrix(FALSE, K, K)), arr.ind = TRUE)
    data.frame(
      k = ij[, 1L],
      l = ij[, 2L],
      key = apply(ij, 1L, function(z) edge_key(z[[1L]], z[[2L]])),
      p_k = p_list[ij[, 1L]],
      p_l = p_list[ij[, 2L]],
      stringsAsFactors = FALSE
    )
  }

  make_incidence_layout <- function(p_list) {
    K <- length(p_list)
    edges <- make_edge_table(p_list)
    lapply(seq_len(K), function(k) {
      neighbors <- setdiff(seq_len(K), k)
      widths <- p_list[neighbors]
      ends <- cumsum(widths)
      starts <- c(1L, head(ends, -1L) + 1L)
      cols <- lapply(seq_along(neighbors), function(j) starts[[j]]:ends[[j]])
      names(cols) <- as.character(neighbors)
      list(
        view = k,
        edge_ids = vapply(neighbors, function(l) {
          which(edges$k == min(k, l) & edges$l == max(k, l))
        }, integer(1L)),
        neighbors = neighbors,
        cols = cols,
        nrow = p_list[[k]],
        ncol = sum(widths)
      )
    })
  }

  empty_edge_list <- function(edge_table, value = 0) {
    out <- lapply(seq_len(nrow(edge_table)), function(e) {
      matrix(value, edge_table$p_k[[e]], edge_table$p_l[[e]])
    })
    names(out) <- edge_table$key
    out
  }

  assemble_Mk <- function(C, k, layout) {
    lk <- layout[[k]]
    pieces <- lapply(seq_along(lk$neighbors), function(j) {
      l <- lk$neighbors[[j]]
      A <- if (!is.null(lk$edge_ids)) C[[lk$edge_ids[[j]]]] else C[[edge_key(k, l)]]
      if (k < l) A else t(A)
    })
    do.call(cbind, pieces)
  }

  assemble_all_M <- function(C, layout) {
    lapply(seq_along(layout), function(k) assemble_Mk(C, k, layout))
  }

  extract_edge_slice <- function(W, k, l, endpoint, layout) {
    stopifnot(k < l, endpoint %in% c(k, l))
    other <- if (endpoint == k) l else k
    cc <- layout[[endpoint]]$cols[[as.character(other)]]
    block <- W[, cc, drop = FALSE]
    if (endpoint == k) block else t(block)
  }

  assemble_full_C <- function(C, p_list) {
    idx <- make_block_indices(p_list)
    edge_table <- make_edge_table(p_list)
    p <- sum(p_list)
    out <- matrix(0, p, p)
    for (e in seq_len(nrow(edge_table))) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      A <- C[[edge_table$key[[e]]]]
      out[idx[[k]], idx[[l]]] <- A
      out[idx[[l]], idx[[k]]] <- t(A)
    }
    out
  }

  split_full_C <- function(C_full, p_list) {
    idx <- make_block_indices(p_list)
    edge_table <- make_edge_table(p_list)
    out <- setNames(vector("list", nrow(edge_table)), edge_table$key)
    for (e in seq_len(nrow(edge_table))) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      out[[edge_table$key[[e]]]] <-
        C_full[idx[[k]], idx[[l]], drop = FALSE]
    }
    out
  }

  center_views <- function(views) {
    means <- lapply(views, colMeans)
    centered <- Map(function(X, m) sweep(X, 2L, m, "-"), views, means)
    list(views = centered, means = means)
  }

  center_views_at <- function(views, means) {
    Map(function(X, m) sweep(X, 2L, m, "-"), views, means)
  }

  # =============================================================================
  # 2. Covariance preparation and objective
  # =============================================================================

  # All fields below depend only on this training/population covariance.
  # A fresh environment is created for every prep object; it never mixes folds.
  cache_problem_matrices <- function(prep) {
    prep$indices <- make_block_indices(prep$p_list)
    prep$Sigma0 <- block_diag(prep$S_kk)
    prep$eig_products <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      outer(prep$eig[[k]]$values, prep$eig[[l]]$values, "*")
    })
    prep$edge_cols_k <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
      k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
      prep$layout[[k]]$cols[[as.character(l)]]
    })
    prep$edge_cols_l <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
      k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
      prep$layout[[l]]$cols[[as.character(k)]]
    })
    prep$loading_factor_cache <- new.env(parent = emptyenv())
    prep
  }

  # Reuse the eigendecomposition for recurring selected supports. Both powers
  # come from ONE decomposition with exactly the previous ridge/eigenvalue floor.
  loading_metric_factors <- function(prep, selected, covariance_ridge) {
    key <- paste0(sprintf("%.17g", covariance_ridge), ":", paste(selected, collapse = ","))
    cache <- prep$loading_factor_cache
    if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    Sigma0 <- prep$Sigma0 %||% block_diag(prep$S_kk)
    S <- Sigma0[selected, selected, drop = FALSE]
    scale_diag <- mean(diag(S))
    if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
    ev <- eigen(symmetrize(S), symmetric = TRUE)
    d <- pmax(ev$values + covariance_ridge * scale_diag, 1e-10)
    ans <- list(
      half = tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors),
      inv_half = tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
    )
    if (is.environment(cache) && length(cache) < LOADING_FACTOR_CACHE_MAX) {
      assign(key, ans, envir = cache)
    }
    ans
  }

  prepare_problem <- function(centered_views) {
    K <- length(centered_views)
    n <- nrow(centered_views[[1L]])
    p_list <- vapply(centered_views, ncol, integer(1L))
    if (any(vapply(centered_views, nrow, integer(1L)) != n)) {
      stop("All views must have the same number of rows.")
    }

    edge_table <- make_edge_table(p_list)
    layout <- make_incidence_layout(p_list)
    S_kk <- lapply(centered_views, function(X) crossprod(X) / n)
    S_kl <- setNames(vector("list", nrow(edge_table)), edge_table$key)

    for (e in seq_len(nrow(edge_table))) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      S_kl[[edge_table$key[[e]]]] <-
        crossprod(centered_views[[k]], centered_views[[l]]) / n
    }

    eig <- lapply(S_kk, function(S) {
      ee <- eigen(symmetrize(S), symmetric = TRUE)
      ee$values <- pmax(ee$values, 0)
      ee
    })

    cache_problem_matrices(list(
      n = n,
      K = K,
      p_list = p_list,
      p = sum(p_list),
      edge_table = edge_table,
      layout = layout,
      S_kk = S_kk,
      S_kl = S_kl,
      eig = eig,
      q = sum(edge_table$p_k * edge_table$p_l)
    ))
  }

  prepare_population_problem <- function(population) {
    p_list <- population$p_list
    edge_table <- make_edge_table(p_list)
    eig <- lapply(population$Sigma_kk, function(S) {
      ee <- eigen(symmetrize(S), symmetric = TRUE)
      ee$values <- pmax(ee$values, 0)
      ee
    })
    cache_problem_matrices(list(
      n = Inf,
      K = length(p_list),
      p_list = p_list,
      p = sum(p_list),
      edge_table = edge_table,
      layout = make_incidence_layout(p_list),
      S_kk = population$Sigma_kk,
      S_kl = population$Sigma_kl,
      eig = eig,
      q = sum(edge_table$p_k * edge_table$p_l)
    ))
  }

  pairwise_loss <- function(C, S_kk, S_ll, S_kl) {
    0.5 * sum(C * (S_kk %*% C %*% S_ll)) - sum(S_kl * C)
  }

  operator_objective <- function(prep, C, rho_e, lambda_g) {
    loss <- 0
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      key <- prep$edge_table$key[[e]]
      loss <- loss + pairwise_loss(
        C[[key]], prep$S_kk[[k]], prep$S_kk[[l]], prep$S_kl[[key]]
      )
    }
    M <- assemble_all_M(C, prep$layout)
    loss + rho_e * sum(vapply(C, function(A) sum(abs(A)), numeric(1L))) +
      lambda_g * sum(vapply(M, function(A) sum(row_l2(A)), numeric(1L)))
  }

  # =============================================================================
  # 3. Oracle1-only zero-penalty consensus solver
  # =============================================================================

  fit_oracle_consensus_admm <- function(
      prep,
      mu_z = 1,
      mu_g = 1,
      max_iter = 2000L,
      abs_tol = 1e-5,
      rel_tol = 1e-4,
      adaptive_mu = TRUE,
      balance_ratio = 10,
      scale_factor = 2,
      adapt_every = 10L,
      group_zero_tol = 1e-8,
      entry_zero_tol = 1e-10,
      init = NULL,
      keep_history = FALSE,
      verbose = FALSE,
      check_every = 1L) {

    # Oracle1 only: both coefficients are fixed at zero, never tuned.
    rho_e <- lambda_g <- 0
    if (mu_z <= 0 || mu_g <= 0) stop("ADMM augmentation parameters must be positive.")

    C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
    Z <- if (!is.null(init$Z)) init$Z else lapply(C, function(A) A)
    H <- if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)

    M0 <- assemble_all_M(C, prep$layout)
    G <- if (!is.null(init$G)) init$G else lapply(M0, function(A) A)
    V <- if (!is.null(init$V)) init$V else lapply(M0, function(A) matrix(0, nrow(A), ncol(A)))

    history <- if (keep_history) {
      data.frame(
        iter = integer(0), objective = numeric(0),
        primal = numeric(0), dual = numeric(0),
        eps_primal = numeric(0), eps_dual = numeric(0),
        mu_z = numeric(0), mu_g = numeric(0)
      )
    } else NULL

    converged <- FALSE
    r_primal <- Inf
    r_dual <- Inf
    eps_primal <- NA_real_
    eps_dual <- NA_real_

    if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
    edge_k <- prep$edge_table$k
    edge_l <- prep$edge_table$l
    edge_keys <- prep$edge_table$key
    n_edges <- length(edge_keys)
    Q <- lapply(prep$eig, `[[`, "vectors")
    denominator_cache <- NULL
    previous_shift <- NA_real_

    for (iter in seq_len(max_iter)) {
      current_shift <- mu_z + 2 * mu_g
      if (!identical(current_shift, previous_shift)) {
        denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
        previous_shift <- current_shift
      }
      Z_old <- Z
      G_old <- G
      C_new <- C  # all edges are overwritten; no zero matrices allocated
      # Each endpoint difference is formed once, not once for every incident edge.
      GV <- Map(function(A, B) A - B, G, V)

      # C updates retain exactly the same Sylvester solution in cached eigenbases.
      for (e in seq_len(n_edges)) {
        k <- edge_k[[e]]; l <- edge_l[[e]]
        target_k <- GV[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        target_l <- t(GV[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        B <- prep$S_kl[[e]] + mu_z * (Z[[e]] - H[[e]]) + mu_g * (target_k + target_l)
        C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
        C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
      }

      M_new <- assemble_all_M(C_new, prep$layout)

      Z_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
        key <- prep$edge_table$key[[e]]
        soft_threshold(C_new[[key]] + H[[key]], rho_e / mu_z)
      }), prep$edge_table$key)

      G_new <- lapply(seq_len(prep$K), function(k) {
        row_group_threshold(M_new[[k]] + V[[k]], lambda_g / mu_g)
      })

      H_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
        key <- prep$edge_table$key[[e]]
        H[[key]] + C_new[[key]] - Z_new[[key]]
      }), prep$edge_table$key)

      V_new <- lapply(seq_len(prep$K), function(k) {
        V[[k]] + M_new[[k]] - G_new[[k]]
      })

      C <- C_new
      Z <- Z_new
      G <- G_new
      H <- H_new
      V <- V_new

      # Residuals, tolerances, adaptive-mu, and convergence are only
      # recomputed every `check_every` iterations (always including the
      # first and last iteration). The ADMM updates above still run every
      # iteration; only this bookkeeping is skipped on non-check iterations,
      # which removes most of the R-level overhead across the very large
      # number of solves performed during cross-validation.
      do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

      if (do_check) {
        r_primal_sq <- 0
        for (key in prep$edge_table$key) {
          r_primal_sq <- r_primal_sq + sum((C_new[[key]] - Z_new[[key]])^2)
        }
        for (k in seq_len(prep$K)) {
          r_primal_sq <- r_primal_sq + sum((M_new[[k]] - G_new[[k]])^2)
        }
        r_primal <- sqrt(r_primal_sq)

        dG <- Map(function(A, B) A - B, G_new, G_old)
        r_dual_sq <- 0
        for (e in seq_len(nrow(prep$edge_table))) {
          k <- prep$edge_table$k[[e]]
          l <- prep$edge_table$l[[e]]
          key <- prep$edge_table$key[[e]]
          dGk <- dG[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
          dGl <- t(dG[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
          S_edge <- mu_z * (Z_new[[key]] - Z_old[[key]]) + mu_g * (dGk + dGl)
          r_dual_sq <- r_dual_sq + sum(S_edge^2)
        }
        r_dual <- sqrt(r_dual_sq)

        norm_Ax_sq <- sum(vapply(C_new, function(A) sum(A^2), numeric(1L))) +
          sum(vapply(M_new, function(A) sum(A^2), numeric(1L)))
        norm_Bz_sq <- sum(vapply(Z_new, function(A) sum(A^2), numeric(1L))) +
          sum(vapply(G_new, function(A) sum(A^2), numeric(1L)))
        d_primal <- prep$q + sum(prep$p_list * (prep$p - prep$p_list))
        eps_primal <- sqrt(d_primal) * abs_tol +
          rel_tol * max(sqrt(norm_Ax_sq), sqrt(norm_Bz_sq))

        dual_adj_sq <- 0
        for (e in seq_len(nrow(prep$edge_table))) {
          k <- prep$edge_table$k[[e]]
          l <- prep$edge_table$l[[e]]
          key <- prep$edge_table$key[[e]]
          Vk <- V_new[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
          Vl <- t(V_new[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
          dual_adj <- mu_z * H_new[[key]] + mu_g * (Vk + Vl)
          dual_adj_sq <- dual_adj_sq + sum(dual_adj^2)
        }
        eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * sqrt(dual_adj_sq)

        if (keep_history) {
          objective_value <- operator_objective(prep, C_new, rho_e, lambda_g)
          history <- rbind(
            history,
            data.frame(
              iter = iter,
              objective = objective_value,
              primal = r_primal,
              dual = r_dual,
              eps_primal = eps_primal,
              eps_dual = eps_dual,
              mu_z = mu_z,
              mu_g = mu_g
            )
          )
        }

        if (!is.finite(r_primal) || !is.finite(r_dual)) {
          warning("ADMM produced a non-finite residual.")
          break
        }

        if (r_primal <= eps_primal && r_dual <= eps_dual) {
          converged <- TRUE
          break
        }

        if (adaptive_mu && iter %% adapt_every == 0L) {
          if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
            mu_z <- mu_z * scale_factor
            mu_g <- mu_g * scale_factor
            H <- lapply(H, function(A) A / scale_factor)
            V <- lapply(V, function(A) A / scale_factor)
          } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
            mu_z <- mu_z / scale_factor
            mu_g <- mu_g / scale_factor
            H <- lapply(H, function(A) A * scale_factor)
            V <- lapply(V, function(A) A * scale_factor)
          }
        }

        if (verbose && (iter == 1L || iter %% 100L == 0L)) {
          cat(sprintf(
            "  iter=%d  primal=%.3e (%.3e)  dual=%.3e (%.3e)\n",
            iter, r_primal, eps_primal, r_dual, eps_dual
          ))
        }
      }
    }

    # Enforce both proximal sparsity patterns in the reported estimate.  This
    # changes only the final tolerance-level disagreement between consensus copies.
    active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
    C_hat <- lapply(Z, function(A) {
      A[abs(A) < entry_zero_tol] <- 0
      A
    })
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      key <- prep$edge_table$key[[e]]
      C_hat[[key]][!active_rows[[k]], ] <- 0
      C_hat[[key]][, !active_rows[[l]]] <- 0
    }

    list(
      C_hat = C_hat,
      C = C,
      Z = Z,
      G = G,
      H = H,
      V = V,
      active_rows = active_rows,
      converged = converged,
      iterations = iter,
      primal_residual = r_primal,
      dual_residual = r_dual,
      eps_primal = eps_primal,
      eps_dual = eps_dual,
      rho_e = rho_e,
      lambda_g = lambda_g,
      mu_z = mu_z,
      mu_g = mu_g,
      history = history
    )
  }

  # =============================================================================
  # 3b. l_21-only consensus ADMM
  # =============================================================================
  # Objective: sum_{k<l} L_kl(C_kl) + lambda_g * sum_k ||M_k(C)||_{2,1}.
  # Only C, G and V are iterated; Z, H and mu_z are absent. The original
  # solver in section 3 is retained unchanged solely for the original
  # Oracle1 path (and the unchanged legacy fit_estimator default).

  fit_l21_admm <- function(
      prep,
      lambda_g,
      mu_g = 1,
      max_iter = 2000L,
      abs_tol = 1e-5,
      rel_tol = 1e-4,
      adaptive_mu = TRUE,
      balance_ratio = 10,
      scale_factor = 2,
      adapt_every = 10L,
      group_zero_tol = 1e-8,
      entry_zero_tol = 1e-10,
      init = NULL,
      keep_history = FALSE,
      verbose = FALSE,
      check_every = 1L) {

    if (length(lambda_g) != 1L || !is.finite(lambda_g) || lambda_g < 0) {
      stop("lambda_g must be a finite, nonnegative scalar.")
    }
    if (length(mu_g) != 1L || !is.finite(mu_g) || mu_g <= 0) {
      stop("mu_g must be a finite, positive scalar.")
    }

    C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)

    M0 <- assemble_all_M(C, prep$layout)
    G <- if (!is.null(init$G)) init$G else lapply(M0, function(A) A)
    V <- if (!is.null(init$V)) init$V else lapply(M0, function(A) matrix(0, nrow(A), ncol(A)))

    history <- if (keep_history) {
      data.frame(
        iter = integer(0), objective = numeric(0),
        primal = numeric(0), dual = numeric(0),
        eps_primal = numeric(0), eps_dual = numeric(0),
        mu_g = numeric(0)
      )
    } else NULL

    converged <- FALSE
    r_primal <- Inf
    r_dual <- Inf
    eps_primal <- NA_real_
    eps_dual <- NA_real_

    if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
    edge_k <- prep$edge_table$k
    edge_l <- prep$edge_table$l
    edge_keys <- prep$edge_table$key
    n_edges <- length(edge_keys)
    Q <- lapply(prep$eig, `[[`, "vectors")
    denominator_cache <- NULL
    previous_shift <- NA_real_

    for (iter in seq_len(max_iter)) {
      # Each edge occurs once at each endpoint: M* M = 2 I.
      current_shift <- 2 * mu_g
      if (!identical(current_shift, previous_shift)) {
        denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
        previous_shift <- current_shift
      }
      G_old <- G
      C_new <- C  # all edges are overwritten; no zero matrices allocated
      # Each endpoint difference is formed once, not once for every incident edge.
      GV <- Map(function(A, B) A - B, G, V)

      # C updates retain exactly the same Sylvester solution in cached eigenbases.
      for (e in seq_len(n_edges)) {
        k <- edge_k[[e]]; l <- edge_l[[e]]
        target_k <- GV[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        target_l <- t(GV[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        B <- prep$S_kl[[e]] + mu_g * (target_k + target_l)
        C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
        C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
      }

      M_new <- assemble_all_M(C_new, prep$layout)

      G_new <- lapply(seq_len(prep$K), function(k) {
        row_group_threshold(M_new[[k]] + V[[k]], lambda_g / mu_g)
      })

      V_new <- lapply(seq_len(prep$K), function(k) {
        V[[k]] + M_new[[k]] - G_new[[k]]
      })

      C <- C_new
      G <- G_new
      V <- V_new

      # Residuals, tolerances, adaptive-mu, and convergence are only
      # recomputed every `check_every` iterations (always including the
      # first and last iteration). The ADMM updates above still run every
      # iteration; only this bookkeeping is skipped on non-check iterations,
      # which removes most of the R-level overhead across the very large
      # number of solves performed during cross-validation.
      do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

      if (do_check) {
        r_primal_sq <- 0
        for (k in seq_len(prep$K)) {
          r_primal_sq <- r_primal_sq + sum((M_new[[k]] - G_new[[k]])^2)
        }
        r_primal <- sqrt(r_primal_sq)

        dG <- Map(function(A, B) A - B, G_new, G_old)
        r_dual_sq <- 0
        for (e in seq_len(nrow(prep$edge_table))) {
          k <- prep$edge_table$k[[e]]
          l <- prep$edge_table$l[[e]]
          key <- prep$edge_table$key[[e]]
          dGk <- dG[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
          dGl <- t(dG[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
          S_edge <- mu_g * (dGk + dGl)
          r_dual_sq <- r_dual_sq + sum(S_edge^2)
        }
        r_dual <- sqrt(r_dual_sq)

        norm_Ax_sq <- sum(vapply(M_new, function(A) sum(A^2), numeric(1L)))
        norm_Bz_sq <- sum(vapply(G_new, function(A) sum(A^2), numeric(1L)))
        d_primal <- sum(prep$p_list * (prep$p - prep$p_list))  # 2 * prep$q
        eps_primal <- sqrt(d_primal) * abs_tol +
          rel_tol * max(sqrt(norm_Ax_sq), sqrt(norm_Bz_sq))

        dual_adj_sq <- 0
        for (e in seq_len(nrow(prep$edge_table))) {
          k <- prep$edge_table$k[[e]]
          l <- prep$edge_table$l[[e]]
          key <- prep$edge_table$key[[e]]
          Vk <- V_new[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
          Vl <- t(V_new[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
          dual_adj <- mu_g * (Vk + Vl)
          dual_adj_sq <- dual_adj_sq + sum(dual_adj^2)
        }
        eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * sqrt(dual_adj_sq)

        if (keep_history) {
          objective_value <- operator_objective(prep, C_new, rho_e = 0, lambda_g = lambda_g)
          history <- rbind(
            history,
            data.frame(
              iter = iter,
              objective = objective_value,
              primal = r_primal,
              dual = r_dual,
              eps_primal = eps_primal,
              eps_dual = eps_dual,
              mu_g = mu_g
            )
          )
        }

        if (!is.finite(r_primal) || !is.finite(r_dual)) {
          warning("ADMM produced a non-finite residual.")
          break
        }

        if (r_primal <= eps_primal && r_dual <= eps_dual) {
          converged <- TRUE
          break
        }

        if (adaptive_mu && iter %% adapt_every == 0L) {
          if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
            mu_g <- mu_g * scale_factor
            V <- lapply(V, function(A) A / scale_factor)
          } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
            mu_g <- mu_g / scale_factor
            V <- lapply(V, function(A) A * scale_factor)
          }
        }

        if (verbose && (iter == 1L || iter %% 100L == 0L)) {
          cat(sprintf(
            "  iter=%d  primal=%.3e (%.3e)  dual=%.3e (%.3e)\n",
            iter, r_primal, eps_primal, r_dual, eps_dual
          ))
        }
      }
    }

    # Return C, masking rows/columns selected out by the endpoint groups.
    # The existing entry_zero_tol is retained ONLY as a numerical cutoff;
    # there is no entrywise penalty, Z variable, or soft-thresholding step.
    active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
    C_hat <- lapply(C, function(A) {
      A[abs(A) < entry_zero_tol] <- 0
      A
    })
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      key <- prep$edge_table$key[[e]]
      C_hat[[key]][!active_rows[[k]], ] <- 0
      C_hat[[key]][, !active_rows[[l]]] <- 0
    }

    list(
      C_hat = C_hat,
      C = C,
      G = G,
      V = V,
      active_rows = active_rows,
      converged = converged,
      iterations = iter,
      primal_residual = r_primal,
      dual_residual = r_dual,
      eps_primal = eps_primal,
      eps_dual = eps_dual,
      rho_e = 0,  # compatibility/reporting only; no entrywise penalty
      lambda_g = lambda_g,
      mu_g = mu_g,
      history = history
    )
  }

  # =============================================================================
  # 4. Faster special case: entrywise-only pairwise EGCAR
  # =============================================================================

  fit_l11_admm <- function(
      prep,
      rho_e,
      mu = 1,
      max_iter = 2000L,
      abs_tol = 1e-5,
      rel_tol = 1e-4,
      adaptive_mu = TRUE,
      balance_ratio = 10,
      scale_factor = 2,
      adapt_every = 10L,
      entry_zero_tol = 1e-10,
      init = NULL,
      verbose = FALSE,
      check_every = 1L) {

    C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
    Z <- if (!is.null(init$Z)) init$Z else empty_edge_list(prep$edge_table)
    H <- if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)
    converged <- FALSE
    r_primal <- Inf
    r_dual <- Inf
    eps_primal <- NA_real_
    eps_dual <- NA_real_

    if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
    edge_k <- prep$edge_table$k
    edge_l <- prep$edge_table$l
    edge_keys <- prep$edge_table$key
    n_edges <- length(edge_keys)
    Q <- lapply(prep$eig, `[[`, "vectors")
    denominator_cache <- NULL
    previous_shift <- NA_real_

    for (iter in seq_len(max_iter)) {
      current_shift <- mu
      if (!identical(current_shift, previous_shift)) {
        denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
        previous_shift <- current_shift
      }
      Z_old <- Z
      C_new <- C  # all edges are overwritten; no zero matrices allocated

      for (e in seq_len(n_edges)) {
        k <- edge_k[[e]]; l <- edge_l[[e]]
        B <- prep$S_kl[[e]] + mu * (Z[[e]] - H[[e]])
        C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
        C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
      }

      Z_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
        key <- prep$edge_table$key[[e]]
        soft_threshold(C_new[[key]] + H[[key]], rho_e / mu)
      }), prep$edge_table$key)

      H_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
        key <- prep$edge_table$key[[e]]
        H[[key]] + C_new[[key]] - Z_new[[key]]
      }), prep$edge_table$key)

      C <- C_new
      Z <- Z_new
      H <- H_new

      # See fit_oracle_consensus_admm for rationale: skip residual/convergence/
      # adaptive-mu bookkeeping on non-check iterations to cut overhead
      # across the very large number of solves used in cross-validation.
      do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

      if (do_check) {
        r_primal <- sqrt(sum(vapply(prep$edge_table$key, function(key) {
          sum((C_new[[key]] - Z_new[[key]])^2)
        }, numeric(1L))))
        r_dual <- mu * sqrt(sum(vapply(prep$edge_table$key, function(key) {
          sum((Z_new[[key]] - Z_old[[key]])^2)
        }, numeric(1L))))

        norm_C <- sqrt(sum(vapply(C_new, function(A) sum(A^2), numeric(1L))))
        norm_Z <- sqrt(sum(vapply(Z_new, function(A) sum(A^2), numeric(1L))))
        eps_primal <- sqrt(prep$q) * abs_tol + rel_tol * max(norm_C, norm_Z)
        norm_dual <- mu * sqrt(sum(vapply(H_new, function(A) sum(A^2), numeric(1L))))
        eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * norm_dual

        if (r_primal <= eps_primal && r_dual <= eps_dual) {
          converged <- TRUE
          break
        }

        if (adaptive_mu && iter %% adapt_every == 0L) {
          if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
            mu <- mu * scale_factor
            H <- lapply(H, function(A) A / scale_factor)
          } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
            mu <- mu / scale_factor
            H <- lapply(H, function(A) A * scale_factor)
          }
        }

        if (verbose && (iter == 1L || iter %% 100L == 0L)) {
          cat(sprintf("  pairwise iter=%d primal=%.3e dual=%.3e\n", iter, r_primal, r_dual))
        }
      }
    }

    C_hat <- lapply(Z, function(A) {
      A[abs(A) < entry_zero_tol] <- 0
      A
    })

    list(
      C_hat = C_hat,
      C = C,
      Z = Z,
      H = H,
      converged = converged,
      iterations = iter,
      primal_residual = r_primal,
      dual_residual = r_dual,
      eps_primal = eps_primal,
      eps_dual = eps_dual,
      rho_e = rho_e,
      lambda_g = 0,
      mu_z = mu
    )
  }

  fit_estimator <- function(prep, rho_e, lambda_g, max_iter, keep_history = FALSE,
                             init = NULL, check_every = 1L, l21_only = FALSE) {
    if (isTRUE(l21_only)) {
      if (length(rho_e) != 1L || !is.finite(rho_e) || rho_e != 0)
        stop("L21-only EGCAR requires rho_e = 0.")
      return(fit_l21_admm(prep, lambda_g, mu_g = MU_G, max_iter = max_iter,
        abs_tol = ABS_TOL, rel_tol = REL_TOL, adaptive_mu = ADAPTIVE_MU,
        group_zero_tol = GROUP_ZERO_TOL, entry_zero_tol = ENTRY_ZERO_TOL,
        init = init, check_every = check_every, keep_history = keep_history))
    }
    if (length(lambda_g) != 1L || !is.finite(lambda_g) || lambda_g != 0)
      stop("L11-only EGCAR requires lambda_g = 0; combined penalties are not supported.")
    fit_l11_admm(prep, rho_e, mu = MU_Z, max_iter = max_iter,
      abs_tol = ABS_TOL, rel_tol = REL_TOL, adaptive_mu = ADAPTIVE_MU,
      entry_zero_tol = ENTRY_ZERO_TOL, init = init, check_every = check_every)
  }

  # =============================================================================
  # 5. Loading extraction, validation score, and metrics
  # =============================================================================

  loading_from_operator <- function(
      prep,
      C,
      rank,
      row_threshold = 1e-4,
      covariance_ridge = 1e-4,
      require_positive = TRUE,
      positive_tol = 1e-10) {

    C_full <- assemble_full_C(C, prep$p_list)
    selected <- which(row_l2(C_full) > row_threshold)
    if (length(selected) < rank) {
      return(list(valid = FALSE, reason = "fewer selected rows than rank"))
    }

    factors <- loading_metric_factors(prep, selected, covariance_ridge)
    S_half <- factors$half
    S_inv_half <- factors$inv_half
    R_sel <- symmetrize(
      S_half %*% C_full[selected, selected, drop = FALSE] %*% S_half
    )
    ee <- eigen(R_sel, symmetric = TRUE)
    if (length(ee$values) < rank) {
      return(list(valid = FALSE, reason = "operator dimension below rank"))
    }
    if (require_positive && ee$values[[rank]] <= positive_tol) {
      return(list(valid = FALSE, reason = "fewer than rank positive eigenvalues"))
    }

    U <- ee$vectors[, seq_len(rank), drop = FALSE]
    L_sel <- S_inv_half %*% U
    L <- matrix(0, prep$p, rank)
    L[selected, ] <- L_sel

    list(
      valid = TRUE,
      L = L,
      U = U,
      selected = selected,
      eigenvalues = ee$values[seq_len(rank)],
      generalized_eigenvalues = 1 + ee$values[seq_len(rank)],
      C_full = C_full
    )
  }

  make_validation_covariance <- function(centered_views) {
    n <- nrow(centered_views[[1L]])
    X <- do.call(cbind, centered_views)
    Sigma <- crossprod(X) / n
    # Diagonal blocks are already present in Sigma; do not traverse X again.
    idx <- make_block_indices(vapply(centered_views, ncol, integer(1L)))
    Sigma0 <- block_diag(lapply(idx, function(ii) Sigma[ii, ii, drop = FALSE]))
    list(Sigma = Sigma, Sigma0 = Sigma0)
  }

  validation_score <- function(L, validation, ridge = 1e-8) {
    if (is.null(L) || any(!is.finite(L))) return(-Inf)
    Q <- symmetrize(crossprod(L, validation$Sigma0 %*% L))
    scale_diag <- mean(diag(Q))
    if (!is.finite(scale_diag) || scale_diag <= 0) return(-Inf)
    ev <- eigen(Q, symmetric = TRUE)
    d <- pmax(ev$values + ridge * scale_diag, 1e-10)
    A <- crossprod(L, validation$Sigma %*% L)
    # tr(Q^(-1/2) A Q^(-1/2)) = sum_j (v_j' A v_j) / d_j.
    # No inverse matrix or two extra matrix-matrix products are constructed.
    score <- sum(colSums(ev$vectors * (A %*% ev$vectors)) / d)
    if (is.finite(score)) score else -Inf
  }

  orthonormal_basis <- function(A, tol = 1e-10) {
    ss <- svd(A, nu = min(dim(A)), nv = 0)
    if (length(ss$d) == 0L || ss$d[[1L]] <= 0) return(NULL)
    keep <- which(ss$d > tol * ss$d[[1L]])
    if (length(keep) == 0L) return(NULL)
    ss$u[, keep, drop = FALSE]
  }

  sine_theta_distance <- function(A, B, rank) {
    QA <- orthonormal_basis(A)
    QB <- orthonormal_basis(B)
    if (is.null(QA) || is.null(QB) || ncol(QA) < rank || ncol(QB) < rank) {
      return(NA_real_)
    }
    QA <- QA[, seq_len(rank), drop = FALSE]
    QB <- QB[, seq_len(rank), drop = FALSE]
    cc <- svd(crossprod(QA, QB), nu = 0, nv = 0)$d
    cc <- pmin(1, pmax(0, cc))
    sqrt(sum(pmax(0, 1 - cc^2)))
  }

  support_metrics <- function(C_full, true_active, threshold = 1e-6) {
    selected <- which(row_l2(C_full) > threshold)
    tp <- length(intersect(selected, true_active))
    fp <- length(setdiff(selected, true_active))
    fn <- length(setdiff(true_active, selected))
    precision <- if (length(selected) == 0L) 0 else tp / length(selected)
    recall <- if (length(true_active) == 0L) 1 else tp / length(true_active)
    fdp <- if (length(selected) == 0L) 0 else fp / length(selected)
    c(precision = precision, recall = recall, fdp = fdp, tp = tp, fp = fp, fn = fn)
  }

  evaluate_method <- function(
      method,
      C,
      loading,
      population,
      rank,
      n,
      rep_id,
      fit_time = NA_real_,
      tune_time = 0,
      rho_e = NA_real_,
      lambda_g = NA_real_,
      c_e = NA_real_,
      c_g = NA_real_,
      converged = NA,
      iterations = NA_integer_,
      status = "ok",
      error_message = NA_character_) {

    C_full <- if (!is.null(C)) assemble_full_C(C, population$p_list) else NULL
    L <- if (!is.null(loading) && isTRUE(loading$valid)) loading$L else NULL

    if (!is.null(C_full)) {
      C_error <- frob(C_full - population$Cstar_full)
      C_relative <- C_error / max(frob(population$Cstar_full), .Machine$double.eps)
      sm <- support_metrics(C_full, population$active_global)
    } else {
      C_error <- C_relative <- NA_real_
      sm <- c(precision = NA, recall = NA, fdp = NA, tp = NA, fp = NA, fn = NA)
    }

    if (!is.null(L)) {
      euclidean_error <- sine_theta_distance(L, population$Lstar, rank)
      sigma0_error <- sine_theta_distance(
        population$Sigma0_half %*% L,
        population$Sigma0_half %*% population$Lstar,
        rank
      )
    } else {
      euclidean_error <- sigma0_error <- NA_real_
    }

    data.frame(
      rep = rep_id,
      rank = rank,
      n = n,
      method = method,
      C_error = C_error,
      C_relative_error = C_relative,
      subspace_euclidean = euclidean_error,
      subspace_sigma0 = sigma0_error,
      support_precision = unname(sm[["precision"]]),
      support_recall = unname(sm[["recall"]]),
      support_fdp = unname(sm[["fdp"]]),
      fit_time = fit_time,
      tune_time = tune_time,
      total_time = fit_time + tune_time,
      rho_e = rho_e,
      lambda_g = lambda_g,
      c_e = c_e,
      c_g = c_g,
      converged = converged,
      iterations = iterations,
      status = status,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }

  # =============================================================================
  # 6. Cross-validation: rho_e for entrywise EGCAR, lambda_g for l21-only EGCAR
  # =============================================================================

  make_folds <- function(n, K = 5L, seed = 1L) {
    set.seed(seed)
    sample(rep(seq_len(K), length.out = n))
  }

  make_fold_objects <- function(views, fold_id) {
    fold_labels <- sort(unique(fold_id))
    lapply(fold_labels, function(f) {
      train_idx <- which(fold_id != f)
      val_idx <- which(fold_id == f)
      train_raw <- lapply(views, function(X) X[train_idx, , drop = FALSE])
      val_raw <- lapply(views, function(X) X[val_idx, , drop = FALSE])
      centered <- center_views(train_raw)
      train_views <- centered$views
      val_views <- center_views_at(val_raw, centered$means)
      list(
        # Retain exactly the training data used to build EGCAR's covariance
        # blocks so external solvers receive the same split and centering.
        fold = f,
        train_views = train_views,
        prep = prepare_problem(train_views),
        validation = make_validation_covariance(val_views)
      )
    })
  }

  # This is the ONLY parallel dispatcher used by the EGCAR and external CV
  # engines. One fold is one task; candidates within that fold remain sequential
  # so the existing caches, warm starts, candidate order, and RNG policy survive.
  parallel_map_candidates <- function(indices, FUN) {
    if (PARALLEL_CV) {
      if (as.integer(future::nbrOfWorkers()) != CV_WORKERS) {
        stop("The shared future plan was changed: all CV methods must use CV_WORKERS workers.")
      }
      future.apply::future_lapply(
        indices, FUN, future.seed = TRUE, future.scheduling = 1
      )
    } else {
      lapply(indices, FUN)
    }
  }

  egcar_positive_cv_grid <- function(x, name = "lambda") {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop(name, " must contain finite, strictly positive coefficients; zero is not a CV candidate.",
           call. = FALSE)
    }
    unique(as.numeric(x))
  }

  cross_validate_penalties <- function(
      full_views, full_prep, fold_objects, rank, rho_e_grid, lambda_g_grid,
      method = c("l11", "l21"), max_iter_cv = 1000L, max_iter_final = 2000L,
      check_every = CHECK_EVERY_ADMM) {
    method <- match.arg(method)
    group <- identical(method, "l21")
    grid <- egcar_positive_cv_grid(if (group) lambda_g_grid else rho_e_grid,
                                    if (group) "lambda_g_grid" else "rho_e_grid")
    nfold <- length(fold_objects)
    if (nfold < 2L) stop("At least two shared folds are required.")
    path_order <- order(grid, decreasing = TRUE)
    tuning_start <- proc.time()[[3L]]
    per_fold <- parallel_map_candidates(seq_len(nfold), function(f) {
      set_blas_threads_one()
      fo <- fold_objects[[f]]
      fo$prep$loading_factor_cache <- new.env(parent = emptyenv())
      previous <- NULL
      rows <- vector("list", length(grid))
      for (j in path_order) {
        notes <- character()
        start <- proc.time()[[3L]]
        one <- tryCatch(withCallingHandlers({
          fit <- fit_estimator(fo$prep, rho_e = if (group) 0 else grid[[j]],
            lambda_g = if (group) grid[[j]] else 0, max_iter = max_iter_cv,
            init = previous, check_every = check_every, l21_only = group)
          loading <- egcar_loading_from_operator(fo$prep, fit$C_hat, rank,
            ROW_THRESHOLD, COVARIANCE_RIDGE, require_positive = TRUE, keep_full_C = FALSE)
          score <- if (isTRUE(loading$valid)) validation_score(loading$L, fo$validation) else -Inf
          list(score = score, converged = fit$converged, iterations = fit$iterations,
            state = fit[if (group) c("C", "G", "V") else c("C", "Z", "H")],
            error = if (!isTRUE(loading$valid)) loading$reason else
              if (!is.finite(score)) "Non-finite validation score." else NA_character_)
        }, warning = function(w) {
          notes <<- unique(c(notes, conditionMessage(w)))
          invokeRestart("muffleWarning")
        }), error = function(e) list(score = -Inf, converged = FALSE,
          iterations = NA_integer_, state = NULL, error = conditionMessage(e)))
        if (!is.null(one$state)) previous <- one$state
        rows[[j]] <- data.frame(candidate = j, lambda = grid[[j]], fold = f,
          n_train = fo$prep$n, score = one$score, loss = -one$score,
          converged = one$converged, iterations = one$iterations,
          elapsed = proc.time()[[3L]] - start, error_message = one$error,
          warning_message = if (length(notes)) paste(notes, collapse = " | ") else NA_character_,
          stringsAsFactors = FALSE)
      }
      do.call(rbind, rows)
    })
    fold_table <- do.call(rbind, per_fold)
    table <- summarize_loading_cv(data.frame(lambda = grid, candidate = seq_along(grid)),
                                   fold_table, nfold)
    table$rho_e <- if (group) 0 else table$lambda
    table$lambda_g <- if (group) table$lambda else 0
    valid <- which(is.finite(table$mean_loss))
    best_index <- if (length(valid)) valid[order(table$mean_loss[valid], -table$lambda[valid])[[1L]]] else NA_integer_
    table$selected <- !is.na(best_index) & seq_len(nrow(table)) == best_index
    tuning_time <- proc.time()[[3L]] - tuning_start
    missing_fit <- list(C_hat = NULL, converged = FALSE, iterations = NA_integer_)
    if (is.na(best_index)) {
      return(list(fit = missing_fit, loading = NULL, cv_table = table, cv_fold_table = fold_table,
        best = NULL, rho_e = if (group) 0 else NA_real_, lambda_g = if (group) NA_real_ else 0,
        fit_time = 0, tuning_time = tuning_time, status = "no_valid_cv",
        error = "No positive CV candidate gave a finite requested-rank score on every fold."))
    }
    best <- table[best_index, , drop = FALSE]
    start <- proc.time()[[3L]]
    fit <- tryCatch(fit_estimator(full_prep, best$rho_e[[1L]], best$lambda_g[[1L]],
      max_iter_final, keep_history = group, l21_only = group),
      error = function(e) list(error = conditionMessage(e)))
    fit_time <- proc.time()[[3L]] - start
    if (is.null(fit$C_hat)) return(list(fit = missing_fit, loading = NULL,
      cv_table = table, cv_fold_table = fold_table, best = best, rho_e = best$rho_e[[1L]],
      lambda_g = best$lambda_g[[1L]], fit_time = fit_time, tuning_time = tuning_time,
      status = "refit_error", error = fit$error))
    loading <- tryCatch(egcar_loading_from_operator(full_prep, fit$C_hat, rank,
      ROW_THRESHOLD, COVARIANCE_RIDGE, require_positive = TRUE),
      error = function(e) list(valid = FALSE, reason = conditionMessage(e)))
    status <- if (!isTRUE(loading$valid)) "invalid_loading" else
      if (!isTRUE(fit$converged)) "not_converged" else "ok"
    list(fit = fit, loading = loading, cv_table = table, cv_fold_table = fold_table,
      best = best, rho_e = best$rho_e[[1L]], lambda_g = best$lambda_g[[1L]],
      tuning_time = tuning_time, fit_time = fit_time, status = status,
      error = if (!isTRUE(loading$valid)) loading$reason else
        if (!isTRUE(fit$converged)) "ADMM did not satisfy its stopping rule." else NA_character_)
  }

  # =============================================================================
  # 7. Positive-semidefinite simulation model
  # =============================================================================

  toeplitz_correlation <- function(p, rho) {
    toeplitz(rho^(0:(p - 1L)))
  }

  metric_normalize <- function(G, Sigma) {
    Gram <- symmetrize(crossprod(G, Sigma %*% G))
    G %*% matrix_power_psd(Gram, -0.5, eig_floor = 1e-12)
  }

  make_population <- function(
      p_list,
      rank,
      active_per_view,
      toeplitz_rho,
      signal,
      seed) {

    set.seed(seed)
    K <- length(p_list)
    if (length(toeplitz_rho) == 1L) toeplitz_rho <- rep(toeplitz_rho, K)
    stopifnot(length(toeplitz_rho) == K, rank <= active_per_view)

    Sigma_kk <- lapply(seq_len(K), function(k) {
      toeplitz_correlation(p_list[[k]], toeplitz_rho[[k]])
    })
    active_local <- lapply(p_list, function(p) seq_len(min(active_per_view, p)))

    U <- lapply(seq_len(K), function(k) {
      G <- matrix(0, p_list[[k]], rank)
      G[active_local[[k]], ] <- matrix(
        rnorm(length(active_local[[k]]) * rank),
        nrow = length(active_local[[k]]),
        ncol = rank
      )
      metric_normalize(G, Sigma_kk[[k]])
    })

    B <- lapply(seq_len(K), function(k) {
      sqrt(signal) * Sigma_kk[[k]] %*% U[[k]]
    })
    Psi <- lapply(seq_len(K), function(k) {
      symmetrize(
        Sigma_kk[[k]] - signal * Sigma_kk[[k]] %*%
          U[[k]] %*% t(U[[k]]) %*% Sigma_kk[[k]]
      )
    })

    edge_table <- make_edge_table(p_list)
    Sigma_kl <- setNames(lapply(seq_len(nrow(edge_table)), function(e) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      B[[k]] %*% t(B[[l]])
    }), edge_table$key)
    Cstar <- setNames(lapply(seq_len(nrow(edge_table)), function(e) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      signal * U[[k]] %*% t(U[[l]])
    }), edge_table$key)

    Sigma0 <- block_diag(Sigma_kk)
    idx <- make_block_indices(p_list)
    Sigma <- Sigma0
    for (e in seq_len(nrow(edge_table))) {
      k <- edge_table$k[[e]]
      l <- edge_table$l[[e]]
      key <- edge_table$key[[e]]
      Sigma[idx[[k]], idx[[l]]] <- Sigma_kl[[key]]
      Sigma[idx[[l]], idx[[k]]] <- t(Sigma_kl[[key]])
    }
    Sigma <- symmetrize(Sigma)

    Sigma0_half <- matrix_power_psd(Sigma0, 0.5)
    Sigma0_inv_half <- matrix_power_psd(Sigma0, -0.5)
    Cstar_full <- assemble_full_C(Cstar, p_list)
    Rstar <- symmetrize(Sigma0_half %*% Cstar_full %*% Sigma0_half)
    ee <- eigen(Rstar, symmetric = TRUE)
    Lstar <- Sigma0_inv_half %*% ee$vectors[, seq_len(rank), drop = FALSE]

    idx <- make_block_indices(p_list)
    active_global <- unlist(lapply(seq_len(K), function(k) idx[[k]][active_local[[k]]]))

    list(
      K = K,
      p_list = p_list,
      p = sum(p_list),
      rank = rank,
      signal = signal,
      Sigma_kk = Sigma_kk,
      Sigma_kl = Sigma_kl,
      Sigma = Sigma,
      Sigma0 = Sigma0,
      Sigma0_half = Sigma0_half,
      U = U,
      B = B,
      Psi = Psi,
      Cstar = Cstar,
      Cstar_full = Cstar_full,
      Lstar = Lstar,
      eigenvalues = ee$values[seq_len(rank)],
      active_local = active_local,
      active_global = active_global
    )
  }

  rmvn_psd <- function(n, Sigma) {
    L <- matrix_power_psd(Sigma, 0.5, eig_floor = 0)
    matrix(rnorm(n * nrow(Sigma)), nrow = n, ncol = nrow(Sigma)) %*% t(L)
  }

  simulate_views <- function(population, n, seed) {
    set.seed(seed)
    Z <- matrix(rnorm(n * population$rank), nrow = n, ncol = population$rank)
    lapply(seq_len(population$K), function(k) {
      Z %*% t(population$B[[k]]) + rmvn_psd(n, population$Psi[[k]])
    })
  }

  # Zero-penalty fit for Oracle1.  This calls the NEW
  # l_11+l_21 consensus ADMM directly, even though rho_e=lambda_g=0, so this
  # oracle baseline uses exactly the same splitting and updates as the proposed
  # joint method, with population covariance blocks in prep.
  fit_zero_penalty_oracle_admm <- function(
      prep,
      rank,
      covariance_ridge,
      max_iter,
      source_label,
      keep_history = TRUE) {

    start <- proc.time()[[3L]]

    out <- tryCatch({
      fit <- fit_oracle_consensus_admm(
        prep = prep,
        mu_z = MU_Z,
        mu_g = MU_G,
        max_iter = max_iter,
        abs_tol = ABS_TOL,
        rel_tol = REL_TOL,
        adaptive_mu = ADAPTIVE_MU,
        group_zero_tol = GROUP_ZERO_TOL,
        entry_zero_tol = ENTRY_ZERO_TOL,
        keep_history = keep_history
      )

      loading <- loading_from_operator(
        prep = prep,
        C = fit$C_hat,
        rank = rank,
        row_threshold = ROW_THRESHOLD,
        covariance_ridge = covariance_ridge,
        require_positive = TRUE
      )

      messages <- character(0)
      status <- "ok"
      if (!isTRUE(fit$converged)) {
        status <- "not_converged"
        messages <- c(messages, sprintf(
          "%s zero-penalty oracle ADMM did not satisfy the stopping rule in %d iterations.",
          source_label, fit$iterations
        ))
      }
      if (!isTRUE(loading$valid)) {
        if (identical(status, "ok")) status <- "invalid_loading"
        messages <- c(messages, paste0(
          "Loading extraction failed: ", loading$reason %||% "unknown reason", "."
        ))
      }

      list(
        C = fit$C_hat,
        loading = loading,
        prep = prep,
        raw_fit = fit,
        status = status,
        error = if (length(messages) == 0L) NA_character_ else paste(messages, collapse = " "),
        converged = fit$converged,
        iterations = fit$iterations,
        rho_e = 0,
        lambda_g = 0
      )
    }, error = function(e) {
      list(
        C = NULL,
        loading = NULL,
        prep = prep,
        raw_fit = NULL,
        status = "failed",
        error = conditionMessage(e),
        converged = FALSE,
        iterations = NA_integer_,
        rho_e = 0,
        lambda_g = 0
      )
    })

    out$time <- proc.time()[[3L]] - start
    out
  }

  # Oracle1: exact population Sigma_kk and Sigma_kl, with rho_e=lambda_g=0.
  # There is no penalty search or cross-validation.  Since neither its inputs nor
  # its coefficients depend on n, it is computed once per replication and rank.
  fit_oracle_population <- function(
      population,
      rank,
      prep = NULL,
      keep_history = TRUE) {

    if (is.null(prep)) prep <- prepare_population_problem(population)
    fit_zero_penalty_oracle_admm(
      prep = prep,
      rank = rank,
      covariance_ridge = 0,
      max_iter = ORACLE1_MAX_ITER,
      source_label = "Population",
      keep_history = keep_history
    )
  }

  # Oracle2: use sample covariance matrices but restrict every block to the
  # known true row supports.
  fit_oracle_support <- function(prep, active_local, ridge = 1e-8) {
    C <- empty_edge_list(prep$edge_table)
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      key <- prep$edge_table$key[[e]]
      sk <- active_local[[k]]
      sl <- active_local[[l]]
      A <- prep$S_kk[[k]][sk, sk, drop = FALSE] + ridge * diag(length(sk))
      B <- prep$S_kk[[l]][sl, sl, drop = FALSE] + ridge * diag(length(sl))
      R <- prep$S_kl[[key]][sk, sl, drop = FALSE]
      C_sub <- t(solve(B, t(solve(A, R))))
      C[[key]][sk, sl] <- C_sub
    }
    C
  }

  # Required initializer dependency closure, bundled from:
  # https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R
  # Inspected 2026-09-10. Function bodies below are retained from that source.
  # These four internal functions do not generate folds or perform CV or TGD.
  # The experiment's common-loss CV and penalized TGD remain authoritative.
  # Copyright and permission notice for this bundled code:
  # MIT License
  # 
  # Copyright (c) 2026 Claire Donnat
  # 
  # Permission is hereby granted, free of charge, to any person obtaining a copy
  # of this software and associated documentation files (the "Software"), to deal
  # in the Software without restriction, including without limitation the rights
  # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  # copies of the Software, and to permit persons to whom the Software is
  # furnished to do so, subject to the following conditions:
  # The above copyright notice and this permission notice shall be included in all
  # copies or substantial portions of the Software.
  # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  # SOFTWARE.

  Soft <- function(a,b){
    if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
    sign(a)*pmax(0,abs(a)-b)
  }

  updatePi <- function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
    C <- Pi + 1/tau*A - nu/tau*B%*%Pi%*%B + nu/tau*sqB%*%(H-Gamma)%*%sqB
    D <- rho/tau
    Soft(C,D)
  }
  updateH <- function(sqB,Gamma,nu,Pi,K){
    temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
    temp <- (temp+t(temp))/2
    ev <- eigen(temp, symmetric = TRUE)
    d <- ev$values

    if(sum(pmin(1,pmax(d,0)))<=K){
      dfinal <- pmin(1,pmax(d,0))
      return(ev$vectors%*%diag(dfinal)%*%t(ev$vectors))
    }
    fr <- function(x) sum(pmin(1,pmax(d-x,0)))
    knots <- unique(c((d-1), d))
    knots <- sort(knots, decreasing=TRUE)
    temp2 <- which(sapply(knots, fr) <= K)
    lentemp <- tail(temp2, 1)
    a <- knots[lentemp]
    b <- knots[lentemp+1]
    fa <- sum(pmin(pmax(d-a,0),1))
    fb <- sum(pmin(pmax(d-b,0),1))
    theta <- a + (b-a) * (K-fa)/(fb-fa)
    dfinal <- pmin(1,pmax(d-theta,0))
    ev$vectors%*%diag(dfinal)%*%t(ev$vectors)
  }
  sgca_init_fixed <- function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    A <- (A + t(A))/2
    B <- (B + t(B))/2
    p <- nrow(B)

    evB <- eigen(B, symmetric = TRUE)
    vals <- pmax(evB$values, 0)
    sqB  <- evB$vectors %*% diag(sqrt(vals), p, p) %*% t(evB$vectors)

    tau <- 4 * nu * (max(vals)^2)
    if (!is.finite(tau) || tau <= 0) tau <- 1

    criteria <- Inf
    iter <- 0
    H <- Pi <- oldPi <- diag(1, p)
    Gamma <- matrix(0, p, p)

    while(criteria > epsilon && iter < maxiter){
      for (j in 1:20){
        Pi <- updatePi(B, sqB, A, H, Gamma, nu, rho, Pi, tau)
      }
      H <- updateH(sqB, Gamma, nu, Pi, K)
      Gamma <- Gamma + (sqB %*% Pi %*% sqB - H) * nu
      criteria <- sqrt(sum((Pi - oldPi)^2))
      oldPi <- Pi
      iter <- iter + 1
      if (trace) cat("iter:", iter, "crit:", criteria, "\n")
    }
    list(Pi=Pi,H=H,Gamma=Gamma,iteration=iter,convergence=criteria)
  }

  # =============================================================================
  # 8. Common-loss CV for SGCA, RGCCA, SGCCA, and MultiCCA
  # =============================================================================
  # All four wrappers minimize exactly -validation_score(L, fo$validation).
  # They consume the same fold_objects as EGCAR. No CV routine generates folds.
  # Each L has p rows and rank columns, in the ORIGINAL variable coordinates.
  # Training and validation variables are centered at TRAINING means, and are
  # not standardized. Package-level variable/block scaling is explicitly off.
  #
  # Sources checked 2026-09-08:
  #   https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R
  #   https://github.com/TannisthaM/SGCA/blob/main/R/tenenhaus_cv.R
  #   https://github.com/TannisthaM/SGCA/blob/main/R/multicca_cv.R
  #   https://rgcca-factory.github.io/RGCCA/reference/rgcca.html
  #   https://github.com/cran/PMA/blob/master/R/MultiCCA.R
  #   Gao & Ma (2023), JMLR 24(135), Algorithm 1:
  #   https://www.jmlr.org/papers/volume24/21-0745/21-0745.pdf
  #
  # SGCA initializer updates are retained; a matrix-cached implementation is the
  # default (FAST_SGCA_INITIALIZER=FALSE restores the original package call).
  # SGCA's TGD below restores the penalized update in Gao--Ma Algorithm 1.
  # Unlike sgca_tgd_safe(), it does NOT B-normalize every iterate, which would
  # annihilate the lambda term. Backtracking is a numerical safeguard added here.
  # RGCCA/SGCCA retain the factorial training scheme; only their CV loss changes.
  # Their block/component signs are resolved on TRAINING observations only.

  validation_loss <- function(L, validation, ridge = 1e-8) {
    -validation_score(L, validation, ridge = ridge)
  }

  validate_loading_matrix <- function(L, p, rank) {
    if (is.null(L)) stop("No loading matrix was returned.")
    L <- as.matrix(L)
    if (!identical(dim(L), c(as.integer(p), as.integer(rank))) ||
        any(!is.finite(L))) {
      stop("The loading matrix must be finite and have dimensions p by rank.")
    }
    ss <- svd(L, nu = 0L, nv = 0L)$d
    if (length(ss) < rank || ss[[1L]] <= 0 ||
        ss[[rank]] <= 1e-10 * ss[[1L]]) {
      stop("The returned loading matrix has numerical rank below the requested rank.")
    }
    L
  }

  # Assemble sample Sigma from EGCAR's already cached covariance blocks.
  # In particular, this uses divisor n, not n-1.
  full_covariance_from_prep <- function(prep) {
    S <- prep$Sigma0 %||% block_diag(prep$S_kk)
    ii <- prep$indices %||% make_block_indices(prep$p_list)
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]
      l <- prep$edge_table$l[[e]]
      R <- prep$S_kl[[prep$edge_table$key[[e]]]]
      S[ii[[k]], ii[[l]]] <- R
      S[ii[[l]], ii[[k]]] <- t(R)
    }
    symmetrize(S)
  }

  named_training_blocks <- function(views) {
    out <- lapply(seq_along(views), function(k) {
      X <- as.matrix(views[[k]])
      storage.mode(X) <- "double"
      if (nrow(X) < 2L || ncol(X) < 1L || any(!is.finite(X))) {
        stop("Training blocks must be finite numeric matrices with at least two rows.")
      }
      colnames(X) <- paste0("block", k, "_V", seq_len(ncol(X)))
      rownames(X) <- paste0("sample", seq_len(nrow(X)))
      X
    })
    names(out) <- paste0("block", seq_along(out))
    out
  }

  # A factorial fit is invariant to independent block/component sign changes;
  # the common signed trace score is not. Fix this indeterminacy BEFORE scoring.
  # For each component, maximize its TRAINING sum of pairwise covariances over
  # block signs, fixing the first sign at +1. This is exact for K <= 12; for
  # larger K use deterministic coordinate ascent on the same training criterion.
  # No held-out observation or population loading is used to choose signs.
  # This is a sign convention, not a rotation/realignment to the true subspace.
  # Assemble training Gram blocks from already computed covariances. scale=1
  # gives X_i'X_j/n; scale=n gives the UNNORMALIZED PMA cross-products.
  training_gram_blocks <- function(prep, scale = 1) {
    G <- lapply(seq_len(prep$K), function(k) vector("list", prep$K))
    for (k in seq_len(prep$K)) G[[k]][[k]] <- scale * prep$S_kk[[k]]
    for (e in seq_len(nrow(prep$edge_table))) {
      k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
      G[[k]][[l]] <- scale * prep$S_kl[[e]]
      G[[l]][[k]] <- t(G[[k]][[l]])
    }
    G
  }

  synchronize_block_signs <- function(weights, training_views, gram_blocks = NULL) {
    if (!ALIGN_EXTERNAL_BLOCK_SIGNS) return(weights)
    K <- length(weights)
    rank <- ncol(weights[[1L]])
    n <- nrow(training_views[[1L]])
    if (K < 2L) return(weights)
    patterns <- NULL
    if (K <= 12L) {
      patterns <- cbind(1, as.matrix(expand.grid(
        rep(list(c(1, -1)), K - 1L), KEEP.OUT.ATTRS = FALSE
      )))
    }
    for (h in seq_len(rank)) {
      if (is.null(gram_blocks)) {
        # Reference/fallback for callers without cached training covariances.
        Y <- do.call(cbind, lapply(seq_len(K), function(k) {
          training_views[[k]] %*% weights[[k]][, h, drop = FALSE]
        }))
        G <- crossprod(Y) / n
      } else {
        G <- matrix(0, K, K)
        for (k in seq_len(K)) {
          wk <- weights[[k]][, h, drop = FALSE]
          for (l in seq.int(k, K)) {
            wl <- weights[[l]][, h, drop = FALSE]
            G[k, l] <- as.numeric(crossprod(wk, gram_blocks[[k]][[l]] %*% wl))
            G[l, k] <- G[k, l]
          }
        }
      }
      if (any(!is.finite(G))) stop("Non-finite training component covariance.")
      if (!is.null(patterns)) {
        vals <- rowSums((patterns %*% G) * patterns)
        signs <- patterns[which.max(vals), ]
      } else {
        signs <- rep(1, K)
        diag(G) <- 0
        for (sweep_id in seq_len(100L)) {
          old <- signs
          for (k in 2:K) {
            z <- sum(G[k, ] * signs)
            if (z != 0) signs[[k]] <- sign(z)
          }
          if (identical(old, signs)) break
        }
      }
      for (k in seq_len(K)) weights[[k]][, h] <- weights[[k]][, h] * signs[[k]]
    }
    weights
  }

  stack_block_weights <- function(weights, context, rank) {
    if (!is.list(weights) || length(weights) != length(context$blocks)) {
      stop("The package did not return one loading matrix per block.")
    }
    # Match block order by names where possible; otherwise retain package order.
    if (!is.null(names(weights)) &&
        all(names(context$blocks) %in% names(weights))) {
      weights <- weights[names(context$blocks)]
    }
    weights <- lapply(seq_along(weights), function(k) {
      W <- as.matrix(weights[[k]])
      if (nrow(W) != context$prep$p_list[[k]] || ncol(W) < rank ||
          any(!is.finite(W))) {
        stop("A package returned an invalid block loading matrix.")
      }
      W[, seq_len(rank), drop = FALSE]
    })
    weights <- synchronize_block_signs(weights, context$blocks, context$gram_blocks)
    validate_loading_matrix(do.call(rbind, weights), context$prep$p, rank)
  }

  # Base-R row binding with a union of columns. Needed because external tuning
  # parameters are not EGCAR's rho_e/lambda_g and must have their own columns.
  bind_rows_fill <- function(...) {
    dfs <- list(...)
    dfs <- Filter(function(d) is.data.frame(d) && nrow(d) > 0L, dfs)
    if (length(dfs) == 0L) return(data.frame())
    cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
    dfs <- lapply(dfs, function(d) {
      for (nm in setdiff(cols, names(d))) d[[nm]] <- rep(NA, nrow(d))
      d[, cols, drop = FALSE]
    })
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    out
  }

  annotate_cv_table <- function(tab, method, rep_id, rank, n, best_candidate = NA_integer_) {
    if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
    # EGCAR continues to maximize its original score internally. These extra
    # reporting columns express the mathematically equivalent minimization.
    if (!"mean_loss" %in% names(tab)) tab$mean_loss <- -tab$mean_score
    if (!"sd_loss" %in% names(tab)) tab$sd_loss <- tab$sd_score
    tab$selected <- !is.na(best_candidate) & tab$candidate == best_candidate
    tab$rep <- rep_id
    tab$rank <- rank
    tab$n <- n
    tab$method <- method
    tab
  }

  summarize_loading_cv <- function(grid, fold_table, nfold) {
    out <- lapply(seq_len(nrow(grid)), function(j) {
      rows <- fold_table[fold_table$candidate == grid$candidate[[j]], , drop = FALSE]
      valid <- is.finite(rows$loss)
      complete <- nrow(rows) == nfold && length(unique(rows$fold)) == nfold && all(valid)
      cv <- data.frame(
        mean_loss = if (complete) mean(rows$loss) else Inf,
        sd_loss = if (complete && nfold > 1L) stats::sd(rows$loss) else NA_real_,
        valid_folds = sum(valid),
        converged_folds = if (all(is.na(rows$converged))) NA_integer_ else
          sum(rows$converged, na.rm = TRUE),
        mean_iterations = if (all(is.na(rows$iterations))) NA_real_ else
          mean(rows$iterations, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      cv$mean_score <- -cv$mean_loss
      cv$sd_score <- cv$sd_loss
      cbind(grid[j, , drop = FALSE], cv)
    })
    do.call(rbind, out)
  }

  # Generic CV engine used by ALL four external wrappers below.
  # prepare_context() is called once per fold; fit_candidate() is called once
  # for every grid row. Contexts may cache TRAINING-only initializers.
  # A candidate is eligible only if every fold has a finite loss. Approximate
  # finite fits are eligible (as in EGCAR); their convergence flags are recorded.
  # Ties retain the first grid row, so each wrapper orders its grid explicitly.
  # All-candidate failure is returned as a failure, not replaced by a hidden fit.
  cross_validate_loading_grid <- function(
      full_views, full_prep, fold_objects, rank, grid,
      prepare_context, fit_candidate, label, seed = 1L,
      parallel_folds = PARALLEL_CV) {

    if (!is.data.frame(grid) || nrow(grid) == 0L) stop("The tuning grid is empty.")
    if (length(fold_objects) < 2L) stop("At least two shared folds are required.")
    if (rank < 1L || rank > full_prep$p) stop("Invalid target rank.")
    if (any(vapply(fold_objects, function(fo) is.null(fo$train_views), logical(1L)))) {
      stop("Use this script's make_fold_objects(): training views must be retained.")
    }
    grid$candidate <- seq_len(nrow(grid))
    Kfold <- length(fold_objects)
    p <- full_prep$p
    seed_for <- function(f) as.integer((as.double(seed) + 1009 * f) %% 2147483646 + 1)
    tuning_start <- proc.time()[[3L]]

    per_fold_fun <- function(f) {
      set_blas_threads_one()
      fo <- fold_objects[[f]]
      context_error <- NULL
      context <- tryCatch(
        prepare_context(fo$train_views, fo$prep, final = FALSE),
        error = function(e) { context_error <<- conditionMessage(e); NULL }
      )
      rows <- vector("list", nrow(grid))
      for (j in seq_len(nrow(grid))) {
        # Equal random initialization across candidates in a given fold. This
        # does not alter the supplied folds. All default solvers use SVD starts.
        set.seed(seed_for(f))
        candidate_start <- proc.time()[[3L]]
        notes <- character(0)
        one <- tryCatch(withCallingHandlers({
          if (is.null(context)) stop(context_error %||% "Training preparation failed.")
          fit <- fit_candidate(context, grid[j, , drop = FALSE], final = FALSE)
          L <- validate_loading_matrix(fit$L, p, rank)
          loss <- validation_loss(L, fo$validation)
          if (!is.finite(loss)) stop("The common held-out loss is non-finite.")
          list(loss = loss, converged = fit$converged %||% NA,
               iterations = fit$iterations %||% NA_integer_, error = NA_character_)
        }, warning = function(w) {
          notes <<- unique(c(notes, conditionMessage(w)))
          invokeRestart("muffleWarning")
        }), error = function(e) {
          list(loss = Inf, converged = FALSE, iterations = NA_integer_,
               error = conditionMessage(e))
        })
        rows[[j]] <- data.frame(
          candidate = grid$candidate[[j]], fold = f,
          n_train = fo$prep$n, loss = one$loss, score = -one$loss,
          converged = one$converged, iterations = one$iterations,
          elapsed = proc.time()[[3L]] - candidate_start,
          error_message = one$error,
          warning_message = if (length(notes)) paste(notes, collapse = " | ") else NA_character_,
          stringsAsFactors = FALSE
        )
      }
      do.call(rbind, rows)
    }

    per_fold <- if (isTRUE(parallel_folds) && PARALLEL_CV) {
      parallel_map_candidates(seq_len(Kfold), per_fold_fun)
    } else {
      lapply(seq_len(Kfold), per_fold_fun)
    }
    fold_table <- do.call(rbind, per_fold)
    cv_table <- summarize_loading_cv(grid, fold_table, Kfold)
    valid <- which(is.finite(cv_table$mean_loss))
    best_index <- if (length(valid)) valid[which.min(cv_table$mean_loss[valid])] else NA_integer_
    cv_table$selected <- !is.na(best_index) & seq_len(nrow(cv_table)) == best_index
    tuning_time <- proc.time()[[3L]] - tuning_start
    # Add parameter columns to the long diagnostics, preserving fold order.
    extra_cols <- setdiff(names(grid), "candidate")
    if (length(extra_cols)) {
      fold_table <- cbind(fold_table, grid[
        match(fold_table$candidate, grid$candidate), extra_cols, drop = FALSE
      ])
      rownames(fold_table) <- NULL
    }

    if (is.na(best_index)) {
      msgs <- unique(stats::na.omit(fold_table$error_message))
      return(list(
        L = NULL, loading = NULL, fit_full = NULL,
        cv_table = cv_table, cv_fold_table = fold_table, best = NULL,
        fit_time = 0, tuning_time = tuning_time, time = tuning_time,
        status = "no_valid_cv", converged = FALSE, iterations = NA_integer_,
        error = paste0("No ", label, " candidate had a finite loss on every fold.",
          if (length(msgs)) paste0(" ", paste(head(msgs, 3L), collapse = " | ")) else "")
      ))
    }

    best <- cv_table[best_index, , drop = FALSE]
    fit_start <- proc.time()[[3L]]
    final_notes <- character(0)
    final_out <- tryCatch(withCallingHandlers({
      # Full-sample centering is done only for the final refit.
      train_full <- center_views(full_views)$views
      context <- prepare_context(train_full, full_prep, final = TRUE)
      set.seed(seed_for(Kfold + 1L))
      fitted <- fit_candidate(context, best, final = TRUE)
      fitted$L <- validate_loading_matrix(fitted$L, p, rank)
      fitted
    }, warning = function(w) {
      final_notes <<- unique(c(final_notes, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) list(L = NULL, error = conditionMessage(e)))
    fit_time <- proc.time()[[3L]] - fit_start

    bad <- is.null(final_out$L)
    conv <- if (bad) FALSE else final_out$converged %||% NA
    status <- if (bad) "refit_error" else if (identical(conv, FALSE)) "not_converged" else "ok"
    list(
      L = final_out$L,
      loading = if (!bad) list(valid = TRUE, L = final_out$L) else NULL,
      fit_full = final_out$fit %||% NULL,
      cv_table = cv_table, cv_fold_table = fold_table, best = best,
      fit_time = fit_time, tuning_time = tuning_time,
      time = fit_time + tuning_time,
      status = status, converged = conv,
      iterations = final_out$iterations %||% NA_integer_,
      error = if (bad) final_out$error else if (identical(conv, FALSE))
        "Final solver reached its iteration limit or did not satisfy its stopping rule." else NA_character_,
      warnings = final_notes
    )
  }

  # ---------------------------- SGCA fitting ---------------------------------

  sgca_hard_rows <- function(U, k) {
    U <- as.matrix(U)
    if (k < ncol(U) || k > nrow(U)) stop("SGCA requires rank <= k <= p.")
    if (k < nrow(U)) {
      keep <- order(rowSums(U * U), decreasing = TRUE)[seq_len(k)]
      U[setdiff(seq_len(nrow(U)), keep), ] <- 0
    }
    U
  }

  sgca_metric_normalize <- function(U, B, tol = 1e-10) {
    G <- symmetrize(crossprod(U, B %*% U))
    ev <- eigen(G, symmetric = TRUE)
    d <- ev$values
    if (any(!is.finite(d)) || d[[1L]] <= 0 ||
        min(d) <= tol * d[[1L]]) stop("SGCA has a rank-deficient training Gram matrix.")
    tcrossprod(sweep(U %*% ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
  }

  # Penalized Gao--Ma TGD, using W = sqrt(lambda) * V for numerical stability.
  # For lambda > 0 this is algebraically Algorithm 1, with adaptive step size:
  #
  # W_new = HT_k[W - 2*eta*(-A W + B W (W' B W - lambda I))].
  #
  # Optimize h(W) = -tr(W' A W)/2 + ||W' B W-lambda I||_F^2/4.
  # Initial W = U0 (lambda I + U0' A U0)^(1/2), U0' B U0 = I.
  # Only the OUTPUT is normalized again. In particular, do not project W'BW
  # onto lambda I every iteration: that would remove the penalty gradient.
  # lambda=0 is excluded (the original penalized formulation then is unbounded).
  # Reused for ALL lambda candidates at the same (training fold, rho, k).
  sgca_prepare_tgd_start <- function(A, B, init, k) {
    U0 <- sgca_metric_normalize(sgca_hard_rows(init, k), B)
    ev <- eigen(symmetrize(crossprod(U0, A %*% U0)), symmetric = TRUE)
    list(vectors = ev$vectors, values = ev$values,
         U0_vectors = U0 %*% ev$vectors)
  }

  # Fast matrix implementation of the supplied repository's sgca_init_fixed().
  # It keeps its 20 inner Pi steps, updateH(), dual update, and stopping rule.
  # Source: TannisthaM/SGCA, R/gao_cv_functions.R (inspected 2026-09-08).
  # B, sqrt(B), tau, and the scaled fixed matrices are prepared once per fold.
  # Within an outer iteration H and Gamma do not change during the 20 Pi steps,
  # so sqrt(B)(H-Gamma)sqrt(B) is computed ONCE, rather than 20 times.
  # No warm start is introduced between rho values.
  sgca_prepare_initializer <- function(A, B, reference, nu = 1) {
    A <- symmetrize(A)
    B <- symmetrize(B)
    h_update <- tryCatch(get("updateH", envir = environment(reference), inherits = TRUE),
                         error = function(e) NULL)
    if (!is.function(h_update)) return(NULL)  # keep compatibility with other versions
    ev <- eigen(B, symmetric = TRUE)
    d <- pmax(ev$values, 0)
    sqB <- tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors)
    tau <- 4 * nu * max(d)^2
    if (!is.finite(tau) || tau <= 0) tau <- 1
    list(p = nrow(B), B = B, sqB = sqB, tau = tau, nu = nu,
         A_scaled = (1 / tau) * A, B_scale = nu / tau, update_H = h_update)
  }

  sgca_init_cached <- function(prepared, rho, K,
                               epsilon = SGCA_INIT_TOL,
                               maxiter = SGCA_MAX_ITER_INIT, trace = FALSE) {
    if (is.null(prepared)) stop("Missing SGCA initializer preparation.")
    z <- prepared
    p <- z$p
    H <- Pi <- oldPi <- diag(1, p)
    Gamma <- matrix(0, p, p)
    criteria <- Inf
    iter <- 0L
    threshold <- rho / z$tau
    while (criteria > epsilon && iter < maxiter) {
      fixed_H <- z$B_scale * (z$sqB %*% (H - Gamma) %*% z$sqB)
      for (j in seq_len(20L)) {
        Pi <- soft_threshold(Pi + z$A_scaled -
          z$B_scale * (z$B %*% Pi %*% z$B) + fixed_H, threshold)
      }
      H <- z$update_H(z$sqB, Gamma, z$nu, Pi, K)
      Gamma <- Gamma + (z$sqB %*% Pi %*% z$sqB - H) * z$nu
      criteria <- sqrt(sum((Pi - oldPi)^2))
      oldPi <- Pi
      iter <- iter + 1L
      if (trace) cat("iter:", iter, "crit:", criteria, "\n")
    }
    list(Pi = Pi, H = H, Gamma = Gamma, iteration = iter, convergence = criteria)
  }

  sgca_tgd_penalized <- function(
      A, B, init, rank, k, lambda,
      eta = SGCA_ETA, max_iter = SGCA_MAX_ITER_TGD,
      tol = SGCA_TGD_TOL, max_backtrack = 50L,
      prepared_start = NULL, matrices_prepared = FALSE) {

    if (!is.finite(lambda) || lambda <= 0) stop("SGCA lambda must be strictly positive.")
    if (!is.finite(eta) || eta <= 0 || max_iter < 1L) stop("Invalid SGCA TGD controls.")
    if (!isTRUE(matrices_prepared)) {
      A <- symmetrize(A)
      B <- symmetrize(B)
    }
    if (is.null(prepared_start)) prepared_start <- sgca_prepare_tgd_start(A, B, init, k)
    st <- prepared_start
    W <- tcrossprod(sweep(st$U0_vectors, 2L,
                           sqrt(pmax(st$values + lambda, 0)), "*"), st$vectors)
    lambda_I <- lambda * diag(rank)
    objective <- function(Z, AZ, BZ) {
      -0.5 * sum(Z * AZ) + 0.25 * sum((crossprod(Z, BZ) - lambda_I)^2)
    }
    AW <- A %*% W
    BW <- B %*% W
    value <- objective(W, AW, BW)
    if (!is.finite(value)) stop("Non-finite SGCA initialization objective.")
    converged <- FALSE
    step <- 2 * eta
    total_backtracks <- 0L
    relative_mapping <- Inf

    for (iter in seq_len(max_iter)) {
      grad <- -AW + BW %*% (crossprod(W, BW) - lambda_I)
      if (any(!is.finite(grad))) stop("Non-finite SGCA gradient.")
      trial_step <- min(2 * eta, 1.25 * step)
      accepted <- FALSE
      for (bt in 0:max_backtrack) {
        Wnew <- sgca_hard_rows(W - trial_step * grad, k)
        if (all(is.finite(Wnew))) {
          AWnew <- A %*% Wnew
          BWnew <- B %*% Wnew
          new_value <- objective(Wnew, AWnew, BWnew)
          diff_sq <- sum((Wnew - W)^2)
          # Sufficient descent for the projected/hard-thresholded step.
          slack <- 1e-12 * max(1, abs(value))
          if (is.finite(new_value) &&
              new_value <= value - diff_sq / (4 * trial_step) + slack) {
            accepted <- TRUE
            break
          }
        }
        trial_step <- trial_step / 2
      }
      if (!accepted) stop("SGCA backtracking could not find a finite descent step.")
      total_backtracks <- total_backtracks + bt
      relative_mapping <- sqrt(diff_sq) / (trial_step * max(1, frob(W)))
      relative_change <- sqrt(diff_sq) / max(1, frob(W))
      W <- Wnew
      AW <- AWnew
      BW <- BWnew
      value <- new_value
      step <- trial_step
      # A small step alone is NOT interpreted as convergence.
      if (relative_mapping <= tol && relative_change <= tol) {
        converged <- TRUE
        break
      }
    }
    L <- sgca_metric_normalize(W, B)
    L <- validate_loading_matrix(L, nrow(A), rank)
    list(L = L, converged = converged, iterations = iter,
         objective_scaled = value, final_step = step,
         relative_gradient_mapping = relative_mapping,
         backtracking_steps = total_backtracks,
         lambda = lambda, k = k)
  }

  # Resolve the bundled initializer; the local CV/TGD routines are unchanged.
  get_sgca_initializer <- function() {
    # Bundled initializer and its dependency closure: no optional namespace lookup.
    sgca_init_fixed
  }

  sgca_common_cv <- function(
      full_views, full_prep, fold_objects, rank,
      k_grid = SGCA_K_GRID, rho_grid = SGCA_RHO_GRID,
      lambda_grid = SGCA_LAMBDA_GRID, seed = 1L,
      parallel_folds = PARALLEL_CV) {

    p <- full_prep$p
    if (any(!is.finite(k_grid)) || any(k_grid != floor(k_grid))) {
      stop("SGCA k_grid must contain finite integers.")
    }
    k_grid <- sort(unique(as.integer(k_grid)))
    k_grid <- k_grid[k_grid >= rank & k_grid <= p]
    if (!length(k_grid)) stop("No SGCA k satisfies rank <= k <= p.")
    if (!length(rho_grid) || any(!is.finite(rho_grid)) || any(rho_grid < 0)) {
      stop("SGCA rho_grid must contain nonnegative, finite DIRECT coefficients.")
    }
    if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
      stop("SGCA lambda_grid must contain positive, finite coefficients.")
    }
    grid <- expand.grid(
      sgca_k = k_grid, sgca_rho = sort(unique(rho_grid), decreasing = TRUE),
      sgca_lambda = sort(unique(lambda_grid)), KEEP.OUT.ATTRS = FALSE
    )
    # Exact CV ties: smaller k, then larger rho, then smaller lambda.
    grid <- grid[order(grid$sgca_k, -grid$sgca_rho, grid$sgca_lambda), , drop = FALSE]
    rownames(grid) <- NULL

    prepare_context <- function(train_views, prep, final) {
      A <- full_covariance_from_prep(prep)
      B <- prep$Sigma0 %||% block_diag(prep$S_kk)
      diag(B) <- diag(B) + SGCA_RIDGE_B
      initializer <- get_sgca_initializer()
      list(
        prep = prep, A = A, B = B, initializer = initializer,
        initializer_prepared = if (FAST_SGCA_INITIALIZER)
          sgca_prepare_initializer(A, B, initializer) else NULL,
        init_cache = new.env(parent = emptyenv()),
        tgd_start_cache = new.env(parent = emptyenv())
      )
    }

    fit_candidate <- function(context, parameters, final) {
      rho <- parameters$sgca_rho[[1L]]
      k <- parameters$sgca_k[[1L]]
      lambda <- parameters$sgca_lambda[[1L]]
      key <- sprintf("rho_%.17g", rho)
      if (!exists(key, envir = context$init_cache, inherits = FALSE)) {
        initialized <- tryCatch({
          z <- if (!is.null(context$initializer_prepared)) {
            sgca_init_cached(context$initializer_prepared, rho = rho, K = rank,
              epsilon = SGCA_INIT_TOL, maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
          } else {
            context$initializer(A = context$A, B = context$B, rho = rho, K = rank,
              nu = 1, epsilon = SGCA_INIT_TOL,
              maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
          }
          Pi <- as.matrix(z$Pi)
          if (any(!is.finite(Pi))) stop("Non-finite SGCA initializer.")
          ss <- svd(Pi, nu = rank, nv = 0L)
          if (length(ss$d) < rank || ss$d[[1L]] <= 0 ||
              ss$d[[rank]] <= 1e-10 * ss$d[[1L]]) {
            stop("The SGCA initializer has numerical rank below the requested rank.")
          }
          U <- sweep(ss$u[, seq_len(rank), drop = FALSE], 2L,
                     sqrt(ss$d[seq_len(rank)]), "*")
          list(U = U, convergence = z$convergence, iteration = z$iteration,
               raw = if (isTRUE(final) && RETAIN_BENCHMARK_FITS) z else NULL,
               error = NULL)
        }, error = function(e) list(U = NULL, convergence = NA_real_, iteration = NA_integer_,
                                    raw = NULL, error = conditionMessage(e)))
        assign(key, initialized, envir = context$init_cache)
      }
      ini <- get(key, envir = context$init_cache, inherits = FALSE)
      if (is.null(ini$U)) stop(ini$error)
      start_key <- paste0(key, "_k", k)
      if (!exists(start_key, envir = context$tgd_start_cache, inherits = FALSE)) {
        st <- tryCatch(list(value = sgca_prepare_tgd_start(context$A, context$B, ini$U, k)),
                       error = function(e) list(error = conditionMessage(e)))
        assign(start_key, st, envir = context$tgd_start_cache)
      }
      st <- get(start_key, envir = context$tgd_start_cache, inherits = FALSE)
      if (is.null(st$value)) stop(st$error)
      tgd <- sgca_tgd_penalized(
        A = context$A, B = context$B, init = ini$U,
        rank = rank, k = k, lambda = lambda,
        prepared_start = st$value, matrices_prepared = TRUE
      )
      init_conv <- if (is.numeric(ini$convergence) && length(ini$convergence) == 1L) {
        is.finite(ini$convergence) && ini$convergence <= SGCA_INIT_TOL
      } else NA
      conv <- if (identical(init_conv, FALSE)) FALSE else tgd$converged
      init_iters <- as.integer(ini$iteration %||% NA_integer_)
      list(
        L = tgd$L,
        fit = if (isTRUE(final) && RETAIN_BENCHMARK_FITS) list(initializer = ini$raw, tgd = tgd,
                   k = k, rho = rho, lambda = lambda,
                   rho_rule = "direct coefficient, unchanged on full-sample refit") else NULL,
        converged = conv, iterations = init_iters + tgd$iterations
      )
    }
    cross_validate_loading_grid(
      full_views, full_prep, fold_objects, rank, grid,
      prepare_context, fit_candidate, label = "SGCA", seed = seed,
      parallel_folds = parallel_folds
    )
  }

  # ------------------------- RGCCA and SGCCA --------------------------------

  rgcca_family_common_cv <- function(
      full_views, full_prep, fold_objects, rank, parameter_grid,
      method = c("rgcca", "sgcca"), seed = 1L,
      parallel_folds = PARALLEL_CV) {

    method <- match.arg(method)
    K <- length(full_prep$p_list)
    if (!requireNamespace("RGCCA", quietly = TRUE)) stop("This benchmark requires RGCCA.")
    if (!length(parameter_grid) || any(!is.finite(parameter_grid))) stop("Invalid RGCCA-family grid.")
    if (method == "rgcca") {
      if (any(parameter_grid < 0 | parameter_grid > 1)) stop("RGCCA tau must be in [0, 1].")
      grid <- data.frame(rgcca_tau = sort(unique(parameter_grid), decreasing = TRUE))
    } else {
      lower <- max(1 / sqrt(full_prep$p_list))
      if (any(parameter_grid < lower - 1e-12 | parameter_grid > 1)) {
        stop("Each tied SGCCA sparsity coefficient must be between max_j(1/sqrt(p_j)) and 1.")
      }
      grid <- data.frame(sgcca_sparsity = sort(unique(pmax(lower, parameter_grid))))
    }
    # Ties: stronger shrinkage (larger tau), or stronger sparsity (smaller bound).
    prepare_context <- function(train_views, prep, final) {
      list(blocks = named_training_blocks(train_views), prep = prep,
           gram_blocks = training_gram_blocks(prep))
    }
    fit_candidate <- function(context, parameters, final) {
      args <- list(
        blocks = context$blocks, connection = 1 - diag(K),
        method = method, ncomp = rep(as.integer(rank), K),
        scheme = RGCCA_SCHEME,
        scale = FALSE, scale_block = FALSE, bias = TRUE,
        init = "svd", verbose = FALSE, quiet = TRUE,
        tol = RGCCA_TOL, n_iter_max = RGCCA_MAX_ITER
      )
      if (method == "rgcca") {
        args$tau <- rep(parameters$rgcca_tau[[1L]], K)
      } else {
        # SGCCA tunes sparsity, NOT tau. This is the actual RGCCA API argument.
        args$sparsity <- rep(parameters$sgcca_sparsity[[1L]], K)
      }
      fit <- do.call(RGCCA::rgcca, args)
      # astar (NOT a) maps the ORIGINAL centered blocks to all deflated scores.
      L <- stack_block_weights(fit$astar, context, rank)
      # Package crit histories provide a limit diagnostic, not a guaranteed
      # optimizer convergence certificate. Do not mark every returned fit TRUE.
      hist <- fit$crit
      lens <- if (is.list(hist)) lengths(hist) else length(hist)
      hit_limit <- length(lens) > 0L && any(lens >= RGCCA_MAX_ITER)
      list(L = L, fit = fit,
           converged = if (hit_limit) FALSE else NA,
           iterations = if (length(lens)) as.integer(sum(lens)) else NA_integer_)
    }
    cross_validate_loading_grid(
      full_views, full_prep, fold_objects, rank, grid,
      prepare_context, fit_candidate, label = toupper(method), seed = seed,
      parallel_folds = parallel_folds
    )
  }

  rgcca_common_cv <- function(
      full_views, full_prep, fold_objects, rank,
      tau_grid = RGCCA_TAU_GRID, seed = 1L, parallel_folds = PARALLEL_CV) {
    rgcca_family_common_cv(
      full_views, full_prep, fold_objects, rank, tau_grid,
      method = "rgcca", seed = seed, parallel_folds = parallel_folds
    )
  }

  sgcca_common_cv <- function(
      full_views, full_prep, fold_objects, rank,
      sparsity_grid = SGCCA_SPARSITY_GRID, seed = 1L, parallel_folds = PARALLEL_CV) {
    rgcca_family_common_cv(
      full_views, full_prep, fold_objects, rank, sparsity_grid,
      method = "sgcca", seed = seed, parallel_folds = parallel_folds
    )
  }

  # ------------------------------- MultiCCA ---------------------------------

  # Matrix-only form of PMA::MultiCCA(type="standard", standardize=FALSE).
  # Source/algorithm: PMA R/MultiCCA.R, Witten & Tibshirani (2009).
  # Uses the INSTALLED PMA BinarySearch, soft, and l2n helpers to retain threshold
  # tolerances and zero-vector handling. Update order, SVD starts, the original
  # (not deflated) stopping objective, the 0.001 stopping threshold, and niter
  # are preserved. Only multiplication order and invariant caching are changed.
  # Let G_ij = X_i' X_j. For a component, previous ws.final are fixed, so
  # D_ij = G_ij - W_i diag(diag(W_i' G_ij W_j)) W_j' is fixed as well.
  # The expensive repeated data operation becomes D_ij %*% w_j.
  multicca_gram_fit <- function(gram, ws_init, penalty, niter = MULTICCA_NITER,
                                ncomponents = 1L, centered_gram = gram) {
    K <- length(gram)
    if (K < 2L || length(ws_init) != K) stop("Invalid MultiCCA Gram blocks.")
    if (length(penalty) == 1L) penalty <- rep(penalty, K)
    if (length(penalty) != K || any(!is.finite(penalty))) stop("Invalid MultiCCA penalties.")
    binary_search <- getFromNamespace("BinarySearch", "PMA")
    soft <- getFromNamespace("soft", "PMA")
    l2n <- getFromNamespace("l2n", "PMA")
    final_w <- lapply(seq_len(K), function(k) {
      if (nrow(ws_init[[k]]) != nrow(gram[[k]][[k]]) ||
          ncol(ws_init[[k]]) < ncomponents) stop("Invalid MultiCCA SVD initializer.")
      matrix(0, nrow(ws_init[[k]]), ncomponents)
    })
    names(final_w) <- names(ws_init)
    cors <- numeric(ncomponents)
    iterations <- integer(ncomponents)
    histories <- vector("list", ncomponents)
    for (h in seq_len(ncomponents)) {
      w <- lapply(ws_init, function(W) W[, h])
      D <- gram
      if (h > 1L) {
        previous <- seq_len(h - 1L)
        for (i in seq_len(K - 1L)) {
          Wi <- final_w[[i]][, previous, drop = FALSE]
          for (j in seq.int(i + 1L, K)) {
            Wj <- final_w[[j]][, previous, drop = FALSE]
            d <- colSums(Wi * (gram[[i]][[j]] %*% Wj))
            D[[i]][[j]] <- gram[[i]][[j]] - tcrossprod(sweep(Wi, 2L, d, "*"), Wj)
            D[[j]][[i]] <- t(D[[i]][[j]])
          }
        }
      }
      curiter <- 1L
      crit_old <- -10
      crit <- -20
      history <- numeric(0)
      # PMA checks this ORIGINAL-data objective BEFORE the next sweep.
      while (curiter <= niter && abs(crit_old - crit) / abs(crit_old) > 0.001 && crit_old != 0) {
        crit_old <- crit
        crit <- 0
        for (i in 2:K) for (j in seq_len(i - 1L)) {
          crit <- crit + as.numeric(crossprod(w[[i]], gram[[i]][[j]] %*% w[[j]]))
        }
        history <- c(history, crit)
        curiter <- curiter + 1L
        for (i in seq_len(K)) {
          total <- 0
          for (j in seq_len(K)[-i]) total <- total + D[[i]][[j]] %*% w[[j]]
          threshold <- binary_search(total, penalty[[i]])
          thresholded <- soft(total, threshold)
          w[[i]] <- thresholded / l2n(thresholded)
        }
      }
      for (k in seq_len(K)) final_w[[k]][, h] <- w[[k]]
      iterations[[h]] <- curiter - 1L
      histories[[h]] <- history
      # Diagnostic only: Pearson correlations, not the common CV loss.
      cor_sum <- 0
      for (i in 2:K) for (j in seq_len(i - 1L)) {
        vi <- as.numeric(crossprod(w[[i]], centered_gram[[i]][[i]] %*% w[[i]]))
        vj <- as.numeric(crossprod(w[[j]], centered_gram[[j]][[j]] %*% w[[j]]))
        co <- if (vi > 0 && vj > 0)
          as.numeric(crossprod(w[[i]], centered_gram[[i]][[j]] %*% w[[j]])) / sqrt(vi * vj) else NA_real_
        if (!is.finite(co)) co <- 0
        cor_sum <- cor_sum + pmin(1, pmax(-1, co))
      }
      cors[[h]] <- cor_sum
    }
    out <- list(ws = final_w, ws.init = ws_init, K = K, call = match.call(),
                type = rep("standard", K), penalty = penalty, cors = cors,
                iterations_per_component = iterations, crit = histories,
                matrix_backend = "cached cross-products")
    class(out) <- "MultiCCA"
    out
  }

  multicca_common_cv <- function(
      full_views, full_prep, fold_objects, rank,
      penalty_grid = MULTICCA_L1_GRID, seed = 1L, parallel_folds = PARALLEL_CV) {

    if (!requireNamespace("PMA", quietly = TRUE)) stop("MultiCCA requires PMA.")
    upper <- min(sqrt(full_prep$p_list))
    if (!length(penalty_grid) || any(!is.finite(penalty_grid)) ||
        any(penalty_grid < 1 | penalty_grid > upper + 1e-12)) {
      stop("PMA MultiCCA's tied L1 bound must lie between 1 and min_j sqrt(p_j).")
    }
    grid <- data.frame(multicca_l1_bound = sort(unique(pmin(upper, penalty_grid))))
    # Ties prefer the smaller L1 bound (more sparsity).
    prepare_context <- function(train_views, prep, final) {
      blocks <- named_training_blocks(train_views)
      ws <- lapply(blocks, function(X) {
        if (rank > min(dim(X))) stop("Too few training dimensions for MultiCCA.")
        # PMA accepts these right singular vectors; cache them across penalties.
        svd(X, nu = 0L, nv = rank)$v[, seq_len(rank), drop = FALSE]
      })
      normalized_gram <- training_gram_blocks(prep)
      raw_gram <- lapply(normalized_gram, function(row) lapply(row, function(G) prep$n * G))
      means <- lapply(blocks, colMeans)
      centered_gram <- lapply(seq_along(blocks), function(i) {
        lapply(seq_along(blocks), function(j) {
          raw_gram[[i]][[j]] - prep$n * tcrossprod(means[[i]], means[[j]])
        })
      })
      list(blocks = blocks, prep = prep, ws_init = ws, gram_blocks = normalized_gram,
           raw_gram = raw_gram, centered_gram = centered_gram)
    }
    fit_candidate <- function(context, parameters, final) {
      backend <- match.arg(MULTICCA_BACKEND, c("gram", "PMA"))
      penalties <- rep(parameters$multicca_l1_bound[[1L]], length(context$blocks))
      fit <- if (backend == "gram") {
        multicca_gram_fit(context$raw_gram, context$ws_init, penalties,
          niter = MULTICCA_NITER, ncomponents = rank, centered_gram = context$centered_gram)
      } else {
        PMA::MultiCCA(xlist = context$blocks, penalty = penalties,
          ws = context$ws_init, type = "standard", ncomponents = rank,
          niter = MULTICCA_NITER, standardize = FALSE, trace = FALSE)
      }
      L <- stack_block_weights(fit$ws, context, rank)
      # Keep previous benchmark diagnostics: an iteration count is not a
      # convergence certificate. Extra Gram-kernel counts are stored in fit.
      list(L = L, fit = fit, converged = NA, iterations = NA_integer_)
    }
    cross_validate_loading_grid(
      full_views, full_prep, fold_objects, rank, grid,
      prepare_context, fit_candidate, label = "MultiCCA", seed = seed,
      parallel_folds = parallel_folds
    )
  }

  # Methods are run sequentially, with the SAME fold-level parallelism as EGCAR.
  # Do not run four methods concurrently and also parallelize their folds.
  # Time = method's complete tuning wall time + full-sample refit wall time.
  run_external_benchmarks <- function(
      views, p_list, rank, seed, fold_objects, full_prep) {

    skipped <- function(label, reason) {
      list(label = label, L = NULL, fit_time = NA_real_, tuning_time = 0,
           time = NA_real_, status = "skipped", error = reason,
           converged = NA, iterations = NA_integer_,
           cv_table = data.frame(), cv_fold_table = data.frame(), best = NULL)
    }
    if (!RUN_EXTERNAL_BENCHMARKS) {
      return(setNames(lapply(EXTERNAL_BENCHMARK_METHODS, skipped,
        reason = "RUN_EXTERNAL_BENCHMARKS is FALSE"), EXTERNAL_BENCHMARK_METHODS))
    }
    package_for <- c(SGCA = NA_character_, RGCCA = "RGCCA", SGCCA = "RGCCA", MultiCCA = "PMA")
    run_one <- function(label) {
      set_blas_threads_one()
      pkg <- unname(package_for[[label]])
      if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
        return(skipped(label, paste("Missing package:", pkg)))
      }
      start <- proc.time()[[3L]]
      cat(sprintf("  %s: common-loss CV on %d shared folds with %d allocated worker(s)...\n",
                  label, length(fold_objects), CV_WORKERS))
      out <- tryCatch({
        f <- switch(label,
          SGCA = sgca_common_cv, RGCCA = rgcca_common_cv,
          SGCCA = sgcca_common_cv, MultiCCA = multicca_common_cv
        )
        f(full_views = views, full_prep = full_prep,
          fold_objects = fold_objects, rank = rank, seed = seed,
          parallel_folds = PARALLEL_CV)
      }, error = function(e) {
        list(L = NULL, fit_time = 0, tuning_time = proc.time()[[3L]] - start,
             status = "error", error = conditionMessage(e),
             converged = FALSE, iterations = NA_integer_,
             cv_table = data.frame(), cv_fold_table = data.frame(), best = NULL)
      })
      out$label <- label
      out$time <- out$fit_time + out$tuning_time
      if (!is.null(out$best)) {
        param_cols <- intersect(names(out$best), c(
          "sgca_k", "sgca_rho", "sgca_lambda", "rgcca_tau",
          "sgcca_sparsity", "multicca_l1_bound"
        ))
        selected_text <- paste(vapply(param_cols, function(nm) {
          paste0(nm, "=", format(out$best[[nm]][[1L]], digits = 5))
        }, character(1L)), collapse = ", ")
        cat(sprintf("    %s; mean loss=%.6g; tune=%.2fs; fit=%.2fs; status=%s\n",
          selected_text, out$best$mean_loss[[1L]], out$tuning_time, out$fit_time, out$status))
      }
      out
    }
    setNames(lapply(EXTERNAL_BENCHMARK_METHODS, run_one), EXTERNAL_BENCHMARK_METHODS)
  }

  # =============================================================================
  # 9. Plotting and checkpointing
  # =============================================================================

  summarize_metric <- function(results, metric, rank_value) {
    d <- results[results$rank == rank_value & is.finite(results[[metric]]), , drop = FALSE]
    if (nrow(d) == 0L) return(NULL)
    groups <- split(d, interaction(d$n, d$method, drop = TRUE))
    out <- do.call(rbind, lapply(groups, function(g) {
      data.frame(
        n = g$n[[1L]],
        method = as.character(g$method[[1L]]),
        mean = mean(g[[metric]], na.rm = TRUE),
        se = if (nrow(g) > 1L) stats::sd(g[[metric]], na.rm = TRUE) / sqrt(nrow(g)) else 0,
        stringsAsFactors = FALSE
      )
    }))
    out$method <- factor(out$method, levels = METHOD_ORDER)
    out
  }

  save_metric_plot <- function(
      results, metric, y_label, rank_value, output_dir,
      log_y = TRUE, title_suffix = "") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
    d <- summarize_metric(results, metric, rank_value)
    if (is.null(d) || nrow(d) == 0L) return(invisible(NULL))
    d$plot_value <- if (log_y) pmax(d$mean, 1e-12) else d$mean

    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = n, y = plot_value, group = method, color = method)
    ) +
      # Only color identifies a method: line type and marker shape are fixed.
      ggplot2::geom_line(linewidth = 0.7, linetype = "solid", na.rm = TRUE) +
      ggplot2::geom_point(size = 2.0, shape = 16, na.rm = TRUE) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_color_manual(values = METHOD_COLORS) +
      ggplot2::labs(
        x = "Sample size n", y = y_label,
        title = paste0(y_label, ", r = ", rank_value, title_suffix),
        color = "Method"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")

    if (log_y) p <- p + ggplot2::scale_y_log10()
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(metric, "_rank_", rank_value, ".pdf")),
      plot = p, width = 9, height = 6, device = "pdf"
    )
    invisible(NULL)
  }

  make_plot_set <- function(results, omit_oracle1 = FALSE) {
    if (omit_oracle1) {
      results <- results[results$method != "Oracle1-population", , drop = FALSE]
      subdir <- "plots_without_oracle1"
      title_suffix <- " (Oracle1 omitted)"
    } else {
      subdir <- "plots_all_methods"
      title_suffix <- " (all methods)"
    }
    output_dir <- file.path(OUT_DIR, subdir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    for (r in sort(unique(results$rank))) {
      # SGCA/RGCCA/SGCCA/MultiCCA estimate loadings rather than C, so they appear in the
      # subspace and time plots; operator/support plots contain methods for which
      # a comparable C estimate is available.
      save_metric_plot(results, "C_relative_error", "Relative operator error",
                       r, output_dir, TRUE, title_suffix)
      save_metric_plot(results, "subspace_euclidean",
                       "Euclidean sine-theta distance",
                       r, output_dir, TRUE, title_suffix)
      save_metric_plot(results, "subspace_sigma0",
                       "Sigma0 sine-theta distance",
                       r, output_dir, TRUE, title_suffix)
      save_metric_plot(results, "total_time", "Total time (seconds)",
                       r, output_dir, TRUE, title_suffix)
      save_metric_plot(results, "support_recall", "Row-support recall",
                       r, output_dir, FALSE, title_suffix)
      save_metric_plot(results, "support_fdp", "Row-support FDP",
                       r, output_dir, FALSE, title_suffix)
    }
    invisible(NULL)
  }

  make_all_plots <- function(results) {
    if (!MAKE_PLOTS || !requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
    make_plot_set(results, omit_oracle1 = FALSE)
    make_plot_set(results, omit_oracle1 = TRUE)
    invisible(NULL)
  }

  # =============================================================================
  # 9b. Post-fit loading visualizations (truth is NEVER used in fitting/CV)
  # =============================================================================
  # Raw loading columns have arbitrary scale, sign, and (for a repeated latent
  # eigenspace) rotation. Compare the global column spaces in original variable
  # coordinates, using Q(L)=L(L'L)^(-1/2). Compute Q stably by a thin SVD.
  # The aligned heatmap uses R=argmin_{R'R=I} ||Qhat R-Qtruth||_F.
  # This is ONE global rotation for the stacked matrix, NOT a separate rotation
  # per view. No variable rows are reordered and no support is imposed.
  # Row importance diag(QQ') and projector differences QQ'-Qtruth Qtruth' are
  # invariant to the original choice of loading basis; they need no alignment.

  loading_euclidean_basis <- function(L, rank) {
    L <- as.matrix(L)
    if (ncol(L) != rank || nrow(L) < rank || any(!is.finite(L))) {
      stop("Non-finite loading matrix or incorrect dimensions.")
    }
    ss <- svd(L, nu = rank, nv = rank)
    if (length(ss$d) < rank || ss$d[[1L]] <= 0 ||
        ss$d[[rank]] <= 1e-10 * ss$d[[1L]]) {
      stop("No valid rank-r loading (numerically rank deficient).")
    }
    # U V' is the polar factor, preserving the raw loading orientation before
    # alignment, unlike using U alone. Its columns are Euclidean orthonormal.
    tcrossprod(ss$u[, seq_len(rank), drop = FALSE],
               ss$v[, seq_len(rank), drop = FALSE])
  }

  plot_loading_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) {
      if (!isTRUE(x$valid) || is.null(x$L)) return(NULL)
      return(as.matrix(x$L))
    }
    as.matrix(x)
  }

  save_multipage_loading_pdf <- function(pages, filename, width = 12, height = 8.5) {
    if (!length(pages)) return(invisible(NULL))
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("loading_plot_", tmpdir = dirname(filename), fileext = ".pdf")
    grDevices::pdf(temporary, width = width, height = height,
                   onefile = TRUE, useDingbats = FALSE)
    device_id <- grDevices::dev.cur()
    on.exit({
      if (device_id %in% grDevices::dev.list()) grDevices::dev.off(device_id)
      if (file.exists(temporary)) unlink(temporary)
    }, add = TRUE)
    for (p in pages) print(p)
    grDevices::dev.off(device_id)
    # Preserve the old complete PDF if plotting failed before this point.
    if (!file.copy(temporary, filename, overwrite = TRUE)) {
      stop("Could not write loading PDF: ", filename)
    }
    invisible(filename)
  }

  make_loading_visualizations <- function(
      loadings, population, rank, n, rep_id, fit_results = data.frame(),
      output_root = file.path(OUT_DIR, "loading_visualizations")) {

    if (!MAKE_PLOTS || !MAKE_LOADING_PLOTS) return(invisible(NULL))
    if (!is.null(LOADING_PLOT_N) && !n %in% LOADING_PLOT_N) return(invisible(NULL))
    if (!is.null(LOADING_PLOT_RANKS) && !rank %in% LOADING_PLOT_RANKS) return(invisible(NULL))
    if (!is.null(LOADING_PLOT_REPS) && !rep_id %in% LOADING_PLOT_REPS) return(invisible(NULL))
    tag <- paste0("rep", rep_id, "_r", rank, "_n", n)
    output_dir <- file.path(output_root, tag)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    methods <- METHOD_ORDER
    p <- population$p
    K <- population$K
    pp <- population$p_list
    idx <- make_block_indices(pp)
    Qtrue <- loading_euclidean_basis(population$Lstar, rank)
    Ptrue <- tcrossprod(Qtrue)
    truth_name <- "Ground truth"
    display_order <- c(truth_name, methods)
    records <- setNames(vector("list", length(display_order)), display_order)
    records[[truth_name]] <- list(
      valid = TRUE, raw = population$Lstar, Q = Qtrue, aligned = Qtrue,
      rotation = diag(rank), projector = Ptrue, projector_difference = matrix(0, p, p),
      importance = rowSums(Qtrue * Qtrue), status = "population truth",
      error = NA_character_, converged = NA, alignment_error = 0, projector_error = 0
    )
    for (label in methods) {
      rr <- if (nrow(fit_results) && "method" %in% names(fit_results))
        fit_results[as.character(fit_results$method) == label, , drop = FALSE] else data.frame()
      status <- if (nrow(rr)) as.character(rr$status[[1L]]) else "not recorded"
      conv <- if (nrow(rr) && "converged" %in% names(rr)) as.logical(rr$converged[[1L]]) else NA
      if (identical(conv, FALSE) && identical(status, "ok")) status <- "not_converged"
      one <- tryCatch({
        L <- plot_loading_matrix(loadings[[label]])
        if (is.null(L)) stop("No valid loading matrix was returned.")
        if (nrow(L) != p) stop("Incorrect number of loading rows.")
        Q <- loading_euclidean_basis(L, rank)
        sv <- svd(crossprod(Q, Qtrue), nu = rank, nv = rank)
        rotation <- tcrossprod(sv$u, sv$v)
        aligned <- Q %*% rotation
        P <- tcrossprod(Q)
        list(valid = TRUE, raw = L, Q = Q, aligned = aligned,
             rotation = rotation, projector = P, projector_difference = P - Ptrue,
             importance = rowSums(Q * Q), status = status, converged = conv,
             alignment_error = frob(aligned - Qtrue),
             projector_error = frob(P - Ptrue), error = NA_character_)
      }, error = function(e) list(
        valid = FALSE, raw = tryCatch(plot_loading_matrix(loadings[[label]]),
                                      error = function(e) NULL),
        Q = matrix(NA_real_, p, rank), aligned = matrix(NA_real_, p, rank),
        rotation = NULL, projector = matrix(NA_real_, p, p),
        projector_difference = matrix(NA_real_, p, p), importance = rep(NA_real_, p),
        status = status, converged = conv, alignment_error = NA_real_,
        projector_error = NA_real_, error = conditionMessage(e)
      ))
      records[[label]] <- one
    }

    summary <- do.call(rbind, lapply(display_order, function(label) {
      z <- records[[label]]
      data.frame(method = label, valid_loading = z$valid, fit_status = z$status,
        converged = z$converged, aligned_basis_error = z$alignment_error,
        projector_error = z$projector_error,
        euclidean_sine_theta_error = z$projector_error / sqrt(2),
        error_message = z$error, stringsAsFactors = FALSE)
    }))
    utils::write.csv(summary, file.path(output_dir, "loading_diagnostics.csv"), row.names = FALSE)
    # Save everything needed to reproduce these displays without refitting.
    saveRDS(list(raw_loadings = loadings, population = population, records = records,
                 rank = rank, n = n, rep_id = rep_id, fit_results = fit_results,
                 normalization = "global Euclidean orthonormal basis",
                 alignment = "one global orthogonal Procrustes rotation for display only"),
            file.path(output_dir, "loading_plot_data.rds"))
    variable <- seq_len(p)
    view <- factor(rep(paste0("View ", seq_len(K)), times = pp),
                   levels = paste0("View ", seq_len(K)))
    local_variable <- unlist(lapply(pp, seq_len), use.names = FALSE)
    truth_importance <- rowSums(Qtrue * Qtrue)
    importance_df <- do.call(rbind, lapply(display_order, function(label) {
      data.frame(method = label, view = view, variable = local_variable,
        global_variable = variable, importance = records[[label]]$importance,
        true_importance = truth_importance,
        true_active = variable %in% population$active_global, stringsAsFactors = FALSE)
    }))
    coefficient_df <- do.call(rbind, lapply(display_order, function(label) {
      raw <- records[[label]]$raw
      if (is.null(raw) || !identical(dim(raw), c(as.integer(p), as.integer(rank)))) {
        raw <- matrix(NA_real_, p, rank)
      }
      data.frame(method = label, view = rep(view, rank),
        variable = rep(local_variable, rank), global_variable = rep(variable, rank),
        component = rep(seq_len(rank), each = p),
        raw = as.vector(raw), normalized = as.vector(records[[label]]$Q),
        aligned = as.vector(records[[label]]$aligned),
        true_normalized = as.vector(Qtrue), stringsAsFactors = FALSE)
    }))
    if (SAVE_LOADING_DATA) {
      utils::write.csv(coefficient_df, file.path(output_dir, "loading_coefficients.csv"), row.names = FALSE)
      utils::write.csv(importance_df, file.path(output_dir, "row_importance.csv"), row.names = FALSE)
    }
    writeLines(c(
      "Loading visualization guide",
      "",
      "Matrices are stacked over views in the original variable order.",
      "Every method and the truth are converted to Q=L(L'L)^(-1/2), by thin SVD.",
      "1. aligned_loading_heatmaps.pdf: compare Qhat R with Qtruth. R is one global",
      "   orthogonal Procrustes rotation. The truth is repeated on every page.",
      "   All pages share a single symmetric loading color scale.",
      "2. row_importance_profiles.pdf: compare diag(Qhat Qhat') with diag(Qtruth Qtruth').",
      "   This is variable participation in the Euclidean subspace, not raw coefficients.",
      "   All methods use solid lines and identical circle markers. Only color differs.",
      "   Pale grey bands mark the known true active variables, for display only.",
      "3. projector_difference_heatmaps.pdf: Qhat Qhat' - Qtruth Qtruth'.",
      "   These differences are invariant to signs, rotations and nonsingular changes",
      "   of loading basis. All pages share a single symmetric difference color scale.",
      "No per-view alignment or population normalization is applied to any fitted model.",
      "Truth is used only after fitting/CV. Nothing in these PDFs changes estimates or scores.",
      "Failures/rank-deficient matrices are displayed as missing, not as zero estimates.",
      "Nonconverged but valid loading matrices are plotted and explicitly labeled.",
      "Raw, normalized, and aligned values are in loading_coefficients.csv and the RDS.",
      "Each output folder represents one simulation (no unaligned averaging across replicates)."
    ), file.path(output_dir, "README.txt"))
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Loading values saved, but ggplot2 is missing; loading PDFs were not generated.")
      return(invisible(summary))
    }

    per_page <- max(1L, as.integer(LOADING_METHODS_PER_PAGE))
    groups <- split(methods, ceiling(seq_along(methods) / per_page))
    palettes <- c(setNames("#111111", truth_name), METHOD_COLORS)
    label_text <- setNames(vapply(display_order, function(label) {
      if (label == truth_name) return(label)
      z <- records[[label]]
      if (!z$valid) return(paste0(label, "\n(no valid loading)"))
      if (!identical(z$status, "ok")) return(paste0(label, "\n", z$status))
      label
    }, character(1L)), display_order)
    plot_subtitle <- paste0("Replication ", rep_id, " | rank ", rank, " | n = ", n)
    common_theme <- ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      legend.position = "bottom"
    )
    tick_variables <- if (max(pp) <= 20L) seq_len(max(pp)) else unique(c(1, pretty(c(1, max(pp)), n = 7)))
    tick_variables <- tick_variables[tick_variables >= 1 & tick_variables <= max(pp)]
    finite_max <- function(x, fallback = 1) {
      x <- x[is.finite(x)]
      if (length(x) && max(abs(x)) > 1e-12) max(abs(x)) else fallback
    }
    loading_limit <- finite_max(coefficient_df$aligned)
    heat_pages <- lapply(groups, function(group) {
      show <- c(truth_name, group)
      d <- coefficient_df[coefficient_df$method %in% show, , drop = FALSE]
      d$method <- factor(d$method, levels = show)
      d$view <- factor(d$view, levels = levels(view))
      invalid <- group[!vapply(records[group], function(z) z$valid, logical(1L))]
      p1 <- ggplot2::ggplot(d, ggplot2::aes(x = component, y = variable, fill = aligned)) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(view ~ method, scales = "free_y", space = "free_y",
          labeller = ggplot2::labeller(method = ggplot2::as_labeller(label_text))) +
        ggplot2::scale_x_continuous(breaks = seq_len(rank)) +
        ggplot2::scale_y_reverse(breaks = tick_variables) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
          midpoint = 0, limits = c(-loading_limit, loading_limit), na.value = "grey90") +
        ggplot2::labs(title = "Estimated loadings aligned to ground truth",
          subtitle = paste0(plot_subtitle, " | normalized global loading bases"),
          x = "Component", y = "Variable within view", fill = "Aligned coefficient",
          caption = "One global orthogonal rotation per method, for display only. Grey panels indicate invalid fits.") +
        common_theme
      if (length(invalid)) {
        missing <- expand.grid(method = invalid, view = levels(view), stringsAsFactors = FALSE)
        missing$component <- (rank + 1) / 2
        missing$variable <- (pp[match(missing$view, levels(view))] + 1) / 2
        missing$method <- factor(missing$method, levels = show)
        p1 <- p1 + ggplot2::geom_text(data = missing,
          ggplot2::aes(x = component, y = variable, label = "No valid\nloading"),
          inherit.aes = FALSE, size = 3)
      }
      p1
    })
    save_multipage_loading_pdf(heat_pages,
      file.path(output_dir, "aligned_loading_heatmaps.pdf"),
      width = max(10, 2.7 * (per_page + 1)), height = max(7.5, 2.7 * K))

    active_bands <- do.call(rbind, lapply(seq_len(K), function(k) {
      a <- population$active_local[[k]]
      data.frame(view = factor(rep(paste0("View ", k), length(a)), levels = levels(view)),
        xmin = a - 0.5, xmax = a + 0.5, ymin = -Inf, ymax = Inf)
    }))
    max_importance <- 1.06 * finite_max(importance_df$importance)
    row_pages <- lapply(groups, function(group) {
      show <- c(truth_name, group)
      d <- importance_df[importance_df$method %in% show, , drop = FALSE]
      d$method <- factor(d$method, levels = show)
      d$view <- factor(d$view, levels = levels(view))
      unavailable <- group[!vapply(records[group], function(z) z$valid, logical(1L))]
      caption <- "Variable participation = diagonal of the Euclidean subspace projector. Grey bands mark true active variables."
      if (length(unavailable)) caption <- paste0(caption, "\nNo valid curve: ", paste(unavailable, collapse = ", "), ".")
      ggplot2::ggplot(d, ggplot2::aes(x = variable, y = importance, color = method, group = method)) +
        ggplot2::geom_rect(data = active_bands,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = "grey93", color = NA, inherit.aes = FALSE) +
        ggplot2::geom_line(linewidth = 0.7, linetype = "solid", na.rm = TRUE) +
        ggplot2::geom_point(size = 1.8, shape = 16, na.rm = TRUE) +
        ggplot2::facet_grid(view ~ ., scales = "free_x") +
        ggplot2::scale_x_continuous(breaks = tick_variables) +
        ggplot2::scale_y_continuous(limits = c(0, max_importance), expand = ggplot2::expansion(mult = c(0, 0.02))) +
        ggplot2::scale_color_manual(values = palettes, breaks = show, labels = label_text[show], drop = FALSE) +
        ggplot2::labs(title = "Which variables contribute to the estimated loading space?",
          subtitle = plot_subtitle, x = "Variable within view", y = "Row importance", color = "Method",
          caption = caption) + common_theme
    })
    save_multipage_loading_pdf(row_pages,
      file.path(output_dir, "row_importance_profiles.pdf"),
      width = 11.5, height = max(7.5, 2.7 * K))

    all_differences <- unlist(lapply(records[methods], `[[`, "projector_difference"), use.names = FALSE)
    difference_limit <- finite_max(all_differences)
    projection_groups <- split(methods, ceiling(seq_along(methods) / 4L))
    boundaries <- head(cumsum(pp), -1L) + 0.5
    view_centers <- cumsum(pp) - (pp - 1) / 2
    projection_pages <- lapply(projection_groups, function(group) {
      d <- do.call(rbind, lapply(group, function(label) {
        ij <- expand.grid(row_variable = seq_len(p), column_variable = seq_len(p), KEEP.OUT.ATTRS = FALSE)
        ij$method <- label
        ij$difference <- as.vector(records[[label]]$projector_difference)
        ij
      }))
      d$method <- factor(d$method, levels = group)
      plot <- ggplot2::ggplot(d, ggplot2::aes(x = column_variable, y = row_variable, fill = difference)) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(~ method, ncol = 2L,
          labeller = ggplot2::labeller(method = ggplot2::as_labeller(label_text))) +
        ggplot2::geom_vline(xintercept = boundaries, color = "grey35", linewidth = 0.3, linetype = "solid") +
        ggplot2::geom_hline(yintercept = boundaries, color = "grey35", linewidth = 0.3, linetype = "solid") +
        ggplot2::scale_x_continuous(breaks = view_centers, labels = paste0("View ", seq_len(K)), expand = c(0, 0)) +
        ggplot2::scale_y_reverse(breaks = view_centers, labels = paste0("View ", seq_len(K)), expand = c(0, 0)) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
          midpoint = 0, limits = c(-difference_limit, difference_limit), na.value = "grey90") +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Estimated minus true subspace projector", subtitle = plot_subtitle,
          x = "Variables grouped by view", y = "Variables grouped by view", fill = "Projector difference",
          caption = "Qhat Qhat' - Qtruth Qtruth'. Invariant to loading scale, sign, rotation and nonsingular basis changes.\nWhite means near-zero discrepancy; grey indicates a missing estimate, not perfect recovery.") + common_theme
      invalid <- group[!vapply(records[group], function(z) z$valid, logical(1L))]
      if (length(invalid)) {
        missing <- data.frame(method = factor(invalid, levels = group),
                              column_variable = (p + 1) / 2, row_variable = (p + 1) / 2)
        plot <- plot + ggplot2::geom_text(data = missing,
          ggplot2::aes(x = column_variable, y = row_variable, label = "No valid loading"),
          inherit.aes = FALSE, size = 3.5)
      }
      plot
    })
    save_multipage_loading_pdf(projection_pages,
      file.path(output_dir, "projector_difference_heatmaps.pdf"), width = 10.5, height = 10.5)
    invisible(summary)
  }

  save_checkpoint <- function(results, cv_results, fits, populations, config,
                              external_cv_fold_results = data.frame()) {
    saveRDS(
      list(
        results = results,
        cv_results = cv_results,
        external_cv_fold_results = external_cv_fold_results,
        fits = fits,
        populations = populations,
        config = config
      ),
      file.path(OUT_DIR, "egcar_simulation.rds")
    )
    utils::write.csv(results, file.path(OUT_DIR, "simulation_results.csv"), row.names = FALSE)
    if (nrow(cv_results) > 0L) {
      utils::write.csv(cv_results, file.path(OUT_DIR, "cv_grid_results.csv"), row.names = FALSE)
      if ("selected" %in% names(cv_results)) {
        chosen <- cv_results[which(cv_results$selected %in% TRUE), , drop = FALSE]
        utils::write.csv(chosen, file.path(OUT_DIR, "selected_cv_parameters.csv"), row.names = FALSE)
      }
    }
    if (SAVE_CV_FOLD_RESULTS && nrow(external_cv_fold_results) > 0L) {
      utils::write.csv(external_cv_fold_results,
        file.path(OUT_DIR, "external_cv_fold_results.csv"), row.names = FALSE)
    }
    if ("status" %in% names(results)) {
      failures <- results[results$status != "ok", , drop = FALSE]
      utils::write.csv(
        failures, file.path(OUT_DIR, "benchmark_failures.csv"), row.names = FALSE
      )
    }
  }


  # Explicit single-penalty fitting entry points for prepared covariance blocks.
  EGCAR_L11 <- function(prep, rho_e, max_iter = MAX_ITER_FINAL, init = NULL,
                         check_every = 1L) {
    fit_estimator(prep, rho_e, 0, max_iter, init = init, check_every = check_every)
  }
  EGCAR_L21 <- function(prep, lambda_g, max_iter = MAX_ITER_FINAL, init = NULL,
                         check_every = 1L, keep_history = FALSE) {
    fit_estimator(prep, 0, lambda_g, max_iter, init = init, check_every = check_every,
                  keep_history = keep_history, l21_only = TRUE)
  }

  # =============================================================================
  # 9c. EGCAR acceleration layer (l11-only and l21-only; no estimator changes)
  # =============================================================================
  # Computational reference: cran/ccar3 R/ecca.r (ecca/ecca_across_lambdas),
  # inspected 2026-09-10, ccar3 DESCRIPTION version 0.1.2.
  # https://github.com/cran/ccar3/blob/master/R/ecca.r
  # Adapted ideas: reduced covariance bases, cached projected cross-covariances,
  # dimension-aware products, reuse along each fold's regularization path.
  # Original extensions here: fused compiled loops, edge-oriented group copies,
  # full-space residuals, spectral objective evaluation, and blockwise GCA loading
  # normalization. ccar3's projected stopping rule, preprocessing, group weights,
  # and two-view loading normalization are NOT substituted for the EGCAR rules.
  #
  # EGCAR_BACKEND=auto       optional RcppArmadillo acceleration, else optimized R
  # EGCAR_BACKEND=cpp        require the compiled backend (fail loudly if unavailable)
  # EGCAR_BACKEND=R          optimized, dependency-free R ADMM backend
  # EGCAR_BACKEND=reference  original l11/l21 solvers and original loading extraction
  # Optional compiled backend: install.packages(c("Rcpp", "RcppArmadillo")); a
  # C++14 compiler is needed. Compilation/loading is timed separately from fitting.
  # There is NO requirement to install ccar3 and NO package installation at runtime.
  #
  # Wide training views use thin SVD bases (all positive singular values retained).
  # Other views reuse the original cached covariance eigenbases, omitting only
  # exact zero eigenvalues. There is NO user-chosen singular-value cutoff. The
  # original covariance and full eig fields remain available to the other methods.
  # Null-space terms are explicitly retained in every Sylvester update.
  # All original grids, folds, penalty formulas, warm-start state conventions,
  # adaptive-mu rules, iteration caps, residual check frequencies and output keys
  # are unchanged. Floating-point reassociation can affect a numerical near-tie.

  EGCAR_BACKEND <- backend
  EGCAR_PARTIAL_EIGEN <- !identical(Sys.getenv("EGCAR_PARTIAL_EIGEN", "1"), "0")
  EGCAR_PARTIAL_EIGEN_MIN <- 128L
  EGCAR_CPP_SOURCE <- paste(c(
    "#include <RcppArmadillo.h>",
    "// [[Rcpp::depends(RcppArmadillo)]]",
    "// [[Rcpp::plugins(cpp14)]]",
    "// Original implementation of the EGCAR ADMM equations. No ccar3 code is copied.",
    "// The spectral reduction and dimension-aware products follow its computational ideas.",
    "#ifndef EGCAR_CORE_HPP",
    "#define EGCAR_CORE_HPP",
    "#include <vector>",
    "#include <array>",
    "#include <cmath>",
    "#include <limits>",
    "#include <stdexcept>",
    "#include <algorithm>",
    "#include <utility>",
    "",
    "namespace egcar_fast {",
    "using arma::mat;",
    "using arma::vec;",
    "struct Edge {",
    "  unsigned k, l;",
    "  mat S, St, D, remainder;",
    "  bool full, project_left, lift_left;",
    "};",
    "struct Problem {",
    "  std::vector<mat> Q;",
    "  std::vector<unsigned> sizes;",
    "  std::vector<Edge> edges;",
    "  double q;",
    "};",
    "struct Control {",
    "  double penalty, mu, abs_tol, rel_tol, balance_ratio, scale_factor;",
    "  int max_iter, check_every, adapt_every;",
    "  bool adaptive, history;",
    "};",
    "struct State {",
    "  std::vector<mat> C, Z, H, Gk, Gl, Vk, Vl;",
    "};",
    "struct Result {",
    "  State state;",
    "  bool converged = false;",
    "  int iterations = 0;",
    "  double primal = std::numeric_limits<double>::infinity();",
    "  double dual = std::numeric_limits<double>::infinity();",
    "  double eps_primal = std::numeric_limits<double>::quiet_NaN();",
    "  double eps_dual = std::numeric_limits<double>::quiet_NaN();",
    "  double mu;",
    "  std::vector<std::array<double,7> > history;",
    "};",
    "inline double sqnorm(const mat& A) { return arma::accu(arma::square(A)); }",
    "inline mat project(const Problem& p, unsigned e, const mat& A) {",
    "  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];",
    "  if (z.project_left) { mat tmp=Qk.t()*A; return tmp*Ql; }",
    "  mat tmp=A*Ql; return Qk.t()*tmp;",
    "}",
    "inline mat lift(const Problem& p, unsigned e, const mat& A) {",
    "  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];",
    "  if (z.lift_left) { mat tmp=Qk*A; return tmp*Ql.t(); }",
    "  mat tmp=A*Ql.t(); return Qk*tmp;",
    "}",
    "// S+shift*T is the right-hand side. The remainder/shift and T terms",
    "// retain the entire null-space component, unlike simply reconstructing Ct.",
    "inline mat c_update(const Problem& p, unsigned e, const mat& T,",
    "                    double shift, const mat& denom, mat& Ct) {",
    "  const Edge& z=p.edges[e];",
    "  mat Pt=project(p,e,T);",
    "  Ct=(z.St+shift*Pt)/denom;",
    "  if (z.full) return lift(p,e,Ct);",
    "  return T+z.remainder/shift+lift(p,e,Ct-Pt);",
    "}",
    "inline double group_sum(const Problem& p, const std::vector<mat>& C) {",
    "  std::vector<vec> norms;",
    "  for (unsigned s:p.sizes) norms.push_back(vec(s,arma::fill::zeros));",
    "  for (unsigned e=0;e<p.edges.size();++e) {",
    "    const Edge& z=p.edges[e];",
    "    norms[z.k]+=arma::sum(arma::square(C[e]),1);",
    "    norms[z.l]+=arma::sum(arma::square(C[e]),0).t();",
    "  }",
    "  double ans=0;",
    "  for (const vec& x:norms) ans+=arma::accu(arma::sqrt(x));",
    "  return ans;",
    "}",
    "inline void threshold_entries(const mat& W, double t, mat& Z) {",
    "  Z.set_size(W.n_rows,W.n_cols);",
    "  for (arma::uword j=0;j<W.n_elem;++j) {",
    "    double x=W[j]; Z[j]=(x>t)?x-t:((x < -t)?x+t:0.0);",
    "  }",
    "}",
    "inline Result solve(const Problem& p, State s, const Control& ctl, bool group,",
    "                    void (*interrupt)()=nullptr,",
    "                    void (*progress)(int,double,double,double,double)=nullptr) {",
    "  if (ctl.max_iter<1 || ctl.check_every<1 || ctl.adapt_every<1 ||",
    "      !(ctl.mu>0) || !std::isfinite(ctl.mu) || ctl.penalty<0 || !std::isfinite(ctl.penalty))",
    "    throw std::invalid_argument(\"Invalid EGCAR solver controls.\");",
    "  const unsigned E=p.edges.size();",
    "  Result out; out.mu=ctl.mu;",
    "  if(ctl.history) out.history.reserve(ctl.max_iter/ctl.check_every+2);",
    "  std::vector<mat> den(E), Ct(E), Wk(E), Wl(E);",
    "  double previous_shift=-1;",
    "  for(int it=1;it<=ctl.max_iter;++it) {",
    "    if(interrupt && (it==1 || it%64==0)) interrupt();",
    "    const bool check=(it==1 || it==ctl.max_iter || it%ctl.check_every==0);",
    "    const double shift=(group?2.0:1.0)*out.mu;",
    "    if(shift!=previous_shift) {",
    "      for(unsigned e=0;e<E;++e) den[e]=p.edges[e].D+shift;",
    "      previous_shift=shift;",
    "    }",
    "    double rp2=0,rd2=0,nc2=0,nz2=0,ny2=0;",
    "    if(!group) {",
    "      for(unsigned e=0;e<E;++e) {",
    "        mat T=s.Z[e]-s.H[e];",
    "        s.C[e]=c_update(p,e,T,shift,den[e],Ct[e]);",
    "        mat W=s.C[e]+s.H[e];",
    "        mat Zn; threshold_entries(W,ctl.penalty/out.mu,Zn);",
    "        mat Hn=W-Zn;",
    "        if(check) {",
    "          rp2+=sqnorm(s.C[e]-Zn); rd2+=sqnorm(Zn-s.Z[e]);",
    "          nc2+=sqnorm(s.C[e]); nz2+=sqnorm(Zn); ny2+=sqnorm(Hn);",
    "        }",
    "        s.Z[e]=std::move(Zn); s.H[e]=std::move(Hn);",
    "      }",
    "      if(check) {",
    "        out.primal=std::sqrt(rp2); out.dual=out.mu*std::sqrt(rd2);",
    "        out.eps_primal=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));",
    "        out.eps_dual=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);",
    "      }",
    "    } else {",
    "      std::vector<vec> normsq;",
    "      for(unsigned size:p.sizes) normsq.push_back(vec(size,arma::fill::zeros));",
    "      for(unsigned e=0;e<E;++e) {",
    "        const Edge& z=p.edges[e];",
    "        mat T=0.5*(s.Gk[e]-s.Vk[e]+s.Gl[e]-s.Vl[e]);",
    "        s.C[e]=c_update(p,e,T,shift,den[e],Ct[e]);",
    "        Wk[e]=s.C[e]+s.Vk[e]; Wl[e]=s.C[e]+s.Vl[e];",
    "        normsq[z.k]+=arma::sum(arma::square(Wk[e]),1);",
    "        normsq[z.l]+=arma::sum(arma::square(Wl[e]),0).t();",
    "      }",
    "      const double threshold=ctl.penalty/out.mu;",
    "      for(vec& x:normsq) for(arma::uword i=0;i<x.n_elem;++i) {",
    "        // Same pmax(norm, .Machine$double.eps) convention as the R original.",
    "        double norm=std::sqrt(x[i]);",
    "        x[i]=threshold<=0?1.0:std::max(0.0,1.0-threshold/std::max(norm,std::numeric_limits<double>::epsilon()));",
    "      }",
    "      for(unsigned e=0;e<E;++e) {",
    "        const Edge& z=p.edges[e];",
    "        mat Gkn=Wk[e]; Gkn.each_col()%=normsq[z.k];",
    "        mat Gln=Wl[e]; Gln.each_row()%=normsq[z.l].t();",
    "        mat Vkn=Wk[e]-Gkn, Vln=Wl[e]-Gln;",
    "        if(check) {",
    "          rp2+=sqnorm(s.C[e]-Gkn)+sqnorm(s.C[e]-Gln);",
    "          rd2+=sqnorm(Gkn-s.Gk[e]+Gln-s.Gl[e]);",
    "          nc2+=2.0*sqnorm(s.C[e]); nz2+=sqnorm(Gkn)+sqnorm(Gln);",
    "          ny2+=sqnorm(Vkn+Vln);",
    "        }",
    "        s.Gk[e]=std::move(Gkn); s.Gl[e]=std::move(Gln);",
    "        s.Vk[e]=std::move(Vkn); s.Vl[e]=std::move(Vln);",
    "      }",
    "      if(check) {",
    "        out.primal=std::sqrt(rp2); out.dual=out.mu*std::sqrt(rd2);",
    "        out.eps_primal=std::sqrt(2.0*p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));",
    "        out.eps_dual=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);",
    "        if(ctl.history) {",
    "          double objective=0;",
    "          for(unsigned e=0;e<E;++e) {",
    "            objective+=0.5*arma::accu(p.edges[e].D%arma::square(Ct[e]))-",
    "                       arma::accu(p.edges[e].S%s.C[e]);",
    "          }",
    "          objective+=ctl.penalty*group_sum(p,s.C);",
    "          out.history.push_back({{double(it),objective,out.primal,out.dual,out.eps_primal,out.eps_dual,out.mu}});",
    "        }",
    "      }",
    "    }",
    "    out.iterations=it;",
    "    if(check) {",
    "      if(!std::isfinite(out.primal)||!std::isfinite(out.dual)) break;",
    "      if(out.primal<=out.eps_primal && out.dual<=out.eps_dual) {",
    "        out.converged=true; break;",
    "      }",
    "      if(ctl.adaptive && it%ctl.adapt_every==0) {",
    "        double factor=1;",
    "        if(out.primal>ctl.balance_ratio*std::max(out.dual,std::numeric_limits<double>::epsilon())) factor=ctl.scale_factor;",
    "        else if(out.dual>ctl.balance_ratio*std::max(out.primal,std::numeric_limits<double>::epsilon())) factor=1.0/ctl.scale_factor;",
    "        if(factor!=1) {",
    "          out.mu*=factor;",
    "          if(group) for(unsigned e=0;e<E;++e) {s.Vk[e]/=factor; s.Vl[e]/=factor;}",
    "          else for(unsigned e=0;e<E;++e) s.H[e]/=factor;",
    "        }",
    "      }",
    "      if(progress && (it==1 || it%100==0))",
    "        progress(it,out.primal,out.dual,out.eps_primal,out.eps_dual);",
    "    }",
    "  }",
    "  out.state=std::move(s); return out;",
    "}",
    "} // namespace egcar_fast",
    "#endif",
    "",
    "",
    "static std::vector<arma::mat> egcar_mat_list(Rcpp::List x) {",
    "  std::vector<arma::mat> ans; ans.reserve(x.size());",
    "  for(int i=0;i<x.size();++i) ans.push_back(Rcpp::as<arma::mat>(x[i]));",
    "  return ans;",
    "}",
    "// Construct R matrices explicitly at the nested-container boundary. Do not",
    "// rely on generic wrapping of std::vector<arma::mat> to retain dim attributes.",
    "// Both Armadillo and R use column-major storage; this is an owning copy, with",
    "// no transpose and no pointer into the soon-to-be-destroyed solver state.",
    "static Rcpp::List egcar_matrix_list_to_R(const std::vector<arma::mat>& values) {",
    "  Rcpp::List out(values.size());",
    "  for (std::size_t i = 0; i < values.size(); ++i) {",
    "    const arma::mat& A = values[i];",
    "    if (A.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max()) ||",
    "        A.n_cols > static_cast<arma::uword>(std::numeric_limits<int>::max()))",
    "      Rcpp::stop(\"Native output matrix dimensions exceed R's matrix limits.\");",
    "    Rcpp::NumericMatrix block(static_cast<int>(A.n_rows),",
    "                              static_cast<int>(A.n_cols));",
    "    if (A.n_elem > 0) std::copy(A.begin(), A.end(), block.begin());",
    "    out[i] = block;",
    "  }",
    "  return out;",
    "}",
    "static void egcar_interrupt() { Rcpp::checkUserInterrupt(); }",
    "static void egcar_progress(int it,double rp,double rd,double ep,double ed) {",
    "  Rcpp::Rcout << \"  iter=\" << it << \" primal=\" << rp << \" (\" << ep << \") dual=\" << rd << \" (\" << ed << \")\\n\";",
    "}",
    "// [[Rcpp::export]]",
    "Rcpp::List egcar_native_solve(Rcpp::List context, Rcpp::List state,",
    "                             Rcpp::List controls, bool group, bool verbose=false) {",
    "  using namespace egcar_fast;",
    "  Problem p; p.Q=egcar_mat_list(context[\"Q\"]);",
    "  Rcpp::IntegerVector sizes=context[\"p_list\"];",
    "  for(int size:sizes) p.sizes.push_back(static_cast<unsigned>(size));",
    "  p.q=Rcpp::as<double>(context[\"q\"]);",
    "  Rcpp::IntegerVector ek=context[\"edge_k\"], el=context[\"edge_l\"];",
    "  Rcpp::List S=context[\"S\"], St=context[\"St\"], D=context[\"D\"], rem=context[\"remainder\"];",
    "  Rcpp::LogicalVector full=context[\"full\"], pl=context[\"project_left\"], ll=context[\"lift_left\"];",
    "  for(int e=0;e<ek.size();++e) {",
    "    Edge z; z.k=ek[e]-1; z.l=el[e]-1;",
    "    z.S=Rcpp::as<mat>(S[e]); z.St=Rcpp::as<mat>(St[e]); z.D=Rcpp::as<mat>(D[e]);",
    "    z.remainder=Rcpp::as<mat>(rem[e]);",
    "    z.full=full[e]; z.project_left=pl[e]; z.lift_left=ll[e];",
    "    p.edges.push_back(std::move(z));",
    "  }",
    "  Control ctl;",
    "  ctl.penalty=Rcpp::as<double>(controls[\"penalty\"]); ctl.mu=Rcpp::as<double>(controls[\"mu\"]);",
    "  ctl.abs_tol=Rcpp::as<double>(controls[\"abs_tol\"]); ctl.rel_tol=Rcpp::as<double>(controls[\"rel_tol\"]);",
    "  ctl.balance_ratio=Rcpp::as<double>(controls[\"balance_ratio\"]); ctl.scale_factor=Rcpp::as<double>(controls[\"scale_factor\"]);",
    "  ctl.max_iter=Rcpp::as<int>(controls[\"max_iter\"]); ctl.check_every=Rcpp::as<int>(controls[\"check_every\"]);",
    "  ctl.adapt_every=Rcpp::as<int>(controls[\"adapt_every\"]); ctl.adaptive=Rcpp::as<bool>(controls[\"adaptive\"]);",
    "  ctl.history=Rcpp::as<bool>(controls[\"history\"]);",
    "  State s; s.C=egcar_mat_list(state[\"C\"]);",
    "  if(group) {",
    "    s.Gk=egcar_mat_list(state[\"Gk\"]); s.Gl=egcar_mat_list(state[\"Gl\"]);",
    "    s.Vk=egcar_mat_list(state[\"Vk\"]); s.Vl=egcar_mat_list(state[\"Vl\"]);",
    "  } else { s.Z=egcar_mat_list(state[\"Z\"]); s.H=egcar_mat_list(state[\"H\"]); }",
    "  Result fit=solve(p,std::move(s),ctl,group,egcar_interrupt,verbose?egcar_progress:nullptr);",
    "  Rcpp::List outstate;",
    "  if(group) outstate=Rcpp::List::create(Rcpp::_ [\"C\"]=egcar_matrix_list_to_R(fit.state.C),Rcpp::_ [\"Gk\"]=egcar_matrix_list_to_R(fit.state.Gk),",
    "      Rcpp::_ [\"Gl\"]=egcar_matrix_list_to_R(fit.state.Gl),Rcpp::_ [\"Vk\"]=egcar_matrix_list_to_R(fit.state.Vk),Rcpp::_ [\"Vl\"]=egcar_matrix_list_to_R(fit.state.Vl));",
    "  else outstate=Rcpp::List::create(Rcpp::_ [\"C\"]=egcar_matrix_list_to_R(fit.state.C),Rcpp::_ [\"Z\"]=egcar_matrix_list_to_R(fit.state.Z),Rcpp::_ [\"H\"]=egcar_matrix_list_to_R(fit.state.H));",
    "  Rcpp::NumericMatrix hist(fit.history.size(),7);",
    "  for(unsigned i=0;i<fit.history.size();++i) for(unsigned j=0;j<7;++j) hist(i,j)=fit.history[i][j];",
    "  if(!std::isfinite(fit.primal)||!std::isfinite(fit.dual)) Rcpp::warning(\"ADMM produced a non-finite residual.\");",
    "  return Rcpp::List::create(Rcpp::_ [\"state\"]=outstate,Rcpp::_ [\"converged\"]=fit.converged,",
    "      Rcpp::_ [\"iterations\"]=fit.iterations,Rcpp::_ [\"primal\"]=fit.primal,Rcpp::_ [\"dual\"]=fit.dual,",
    "      Rcpp::_ [\"eps_primal\"]=fit.eps_primal,Rcpp::_ [\"eps_dual\"]=fit.eps_dual,",
    "      Rcpp::_ [\"mu\"]=fit.mu,Rcpp::_ [\"history\"]=hist, Rcpp::_ [\"matrix_api\"]=2);",
    "}"
  ), collapse = "\n")
  EGCAR_CPP_CACHE <- file.path(normalizePath(OUT_DIR, mustWork = TRUE), ".egcar_native_cache")
  EGCAR_CPP_FILE <- file.path(EGCAR_CPP_CACHE, "egcar_native_v2_matrix_api.cpp")

  # Keep the original implementations callable for checks and exact rollback.
  fit_l11_admm_reference <- fit_l11_admm
  fit_l21_admm_reference <- fit_l21_admm
  cache_problem_matrices_reference <- cache_problem_matrices
  prepare_problem_reference <- prepare_problem

  # A compiled function pointer must never be serialized to a multisession worker.
  # Each process obtains its own locally loaded function from the shared build
  # cache; the pid guard also handles a future plan change or a forked process.
  egcar_native_function <- function(required = FALSE) {
    runtime <- getOption("egcar.native.runtime.v2.matrix.api")
    if (is.list(runtime) && identical(runtime$pid, Sys.getpid()) &&
        identical(runtime$file, EGCAR_CPP_FILE)) {
      if (required && is.null(runtime$fun)) stop(runtime$error)
      return(runtime$fun)
    }
    failure <- NULL
    fun <- tryCatch({
      if (!requireNamespace("Rcpp", quietly = TRUE) ||
          !requireNamespace("RcppArmadillo", quietly = TRUE)) {
        stop("The compiled EGCAR backend requires Rcpp and RcppArmadillo.")
      }
      dir.create(EGCAR_CPP_CACHE, recursive = TRUE, showWarnings = FALSE)
      # Master writes this before dispatching any worker, avoiding write races.
      if (!file.exists(EGCAR_CPP_FILE)) writeLines(EGCAR_CPP_SOURCE, EGCAR_CPP_FILE)
      local <- new.env(parent = baseenv())
      Rcpp::sourceCpp(file = EGCAR_CPP_FILE, env = local,
                     cacheDir = EGCAR_CPP_CACHE, rebuild = FALSE,
                     showOutput = FALSE, verbose = FALSE)
      get("egcar_native_solve", envir = local, inherits = FALSE)
    }, error = function(e) { failure <<- conditionMessage(e); NULL })
    options(egcar.native.runtime.v2.matrix.api = list(pid = Sys.getpid(), file = EGCAR_CPP_FILE,
                                         fun = fun, error = failure))
    if (required && is.null(fun)) stop(failure)
    fun
  }

  egcar_project <- function(z, e, A) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    if (z$project_left[[e]]) crossprod(z$Q[[k]], A) %*% z$Q[[l]] else
      crossprod(z$Q[[k]], A %*% z$Q[[l]])
  }

  egcar_lift <- function(z, e, A) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    if (z$lift_left[[e]]) tcrossprod(z$Q[[k]] %*% A, z$Q[[l]]) else
      z$Q[[k]] %*% tcrossprod(A, z$Q[[l]])
  }

  egcar_prepare_context <- function(prep, centered_views = NULL) {
    spectra <- lapply(seq_len(prep$K), function(k) {
      if (!is.null(centered_views) && nrow(centered_views[[k]]) < prep$p_list[[k]]) {
        # Statistical covariance is still crossprod(X)/n. Thin SVD avoids
        # propagating numerical eigenvalues in the algebraic covariance nullspace.
        X <- centered_views[[k]]
        ss <- svd(X, nu = 0L, nv = min(dim(X)))
        keep <- which(ss$d > 0)
        return(list(vectors = ss$v[, keep, drop = FALSE], values = ss$d[keep]^2 / nrow(X),
                    origin = "thin training-data SVD"))
      }
      ev <- prep$eig[[k]]; keep <- which(ev$values > 0)
      list(vectors = ev$vectors[, keep, drop = FALSE], values = ev$values[keep],
           origin = "cached covariance eigendecomposition")
    })
    Q <- lapply(spectra, `[[`, "vectors")
    d <- lapply(spectra, `[[`, "values")
    rk <- lengths(d)
    ek <- prep$edge_table$k; el <- prep$edge_table$l
    pk <- prep$p_list[ek]; pl <- prep$p_list[el]
    a <- rk[ek]; b <- rk[el]
    # Compare total multiply counts, not just the size of one intermediate.
    z <- list(Q = Q, p_list = as.integer(prep$p_list), edge_k = as.integer(ek),
      edge_l = as.integer(el), keys = prep$edge_table$key, q = prep$q,
      cols_k = prep$edge_cols_k, cols_l = prep$edge_cols_l,
      S = unname(prep$S_kl), full = (a == pk & b == pl),
      basis_origin = vapply(spectra, `[[`, character(1L), "origin"),
      project_left = (a * pk * pl + a * pl * b <= pk * pl * b + a * pk * b),
      lift_left = (pk * a * b + pk * b * pl <= a * b * pl + pk * a * pl))
    z$D <- lapply(seq_along(ek), function(e) outer(d[[ek[[e]]]], d[[el[[e]]]], "*"))
    z$St <- lapply(seq_along(ek), function(e) egcar_project(z, e, z$S[[e]]))
    z$remainder <- lapply(seq_along(ek), function(e) {
      if (z$full[[e]]) matrix(0, 0L, 0L) else z$S[[e]] - egcar_lift(z, e, z$St[[e]])
    })
    z
  }

  # Statistics and original covariance/eigendecomposition fields are untouched.
  # The additional EGCAR-only context is prepared once with the shared problem,
  # outside the fit/CV timers, just like the pre-existing eig_products cache.
  cache_problem_matrices <- function(prep) {
    prep <- cache_problem_matrices_reference(prep)
    # The data-aware prepare_problem wrapper handles wide training views.
    if (EGCAR_BACKEND != "reference" && !any(prep$p_list > prep$n))
      prep$egcar_context <- egcar_prepare_context(prep)
    prep
  }

  prepare_problem <- function(centered_views) {
    prep <- prepare_problem_reference(centered_views)
    if (EGCAR_BACKEND != "reference" && is.null(prep$egcar_context))
      prep$egcar_context <- egcar_prepare_context(prep, centered_views)
    prep
  }

  egcar_get_context <- function(prep) {
    if (!is.null(prep$egcar_context)) return(prep$egcar_context)
    if (is.null(prep$edge_cols_k)) prep <- cache_problem_matrices_reference(prep)
    egcar_prepare_context(prep)
  }

  egcar_initial_state <- function(prep, z, init, group) {
    C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
    # Match the original default initialization of every consensus copy.
    if (!group) return(list(C = C,
      Z = if (!is.null(init$Z)) init$Z else empty_edge_list(prep$edge_table),
      H = if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)))
    G <- if (!is.null(init$G)) init$G else assemble_all_M(C, prep$layout)
    V <- if (!is.null(init$V)) init$V else lapply(G, function(A) matrix(0, nrow(A), ncol(A)))
    E <- seq_along(z$edge_k)
    list(C = C,
      Gk = lapply(E, function(e) G[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
      Gl = lapply(E, function(e) t(G[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])),
      Vk = lapply(E, function(e) V[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
      Vl = lapply(E, function(e) t(V[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])))
  }

  egcar_view_copies <- function(z, left, right) {
    p <- sum(z$p_list)
    out <- lapply(z$p_list, function(pk) matrix(0, pk, p - pk))
    for (e in seq_along(z$edge_k)) {
      out[[z$edge_k[[e]]]][, z$cols_k[[e]]] <- left[[e]]
      out[[z$edge_l[[e]]]][, z$cols_l[[e]]] <- t(right[[e]])
    }
    out
  }

  egcar_group_norms <- function(z, left, right = left) {
    ans <- lapply(z$p_list, numeric)
    for (e in seq_along(z$edge_k)) {
      k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
      ans[[k]] <- ans[[k]] + rowSums(left[[e]] * left[[e]])
      ans[[l]] <- ans[[l]] + colSums(right[[e]] * right[[e]])
    }
    lapply(ans, sqrt)
  }

  # Optimized base-R backend; it is not the old reference solver. The native
  # backend below implements these same updates with fused in-place loops.
  egcar_solve_R <- function(z, s, ctl, group, verbose = FALSE) {
    E <- seq_along(z$edge_k); K <- seq_along(z$p_list)
    mu <- ctl$mu; last_shift <- NA_real_; den <- NULL
    converged <- FALSE; rp <- rd <- Inf; ep <- ed <- NA_real_
    history <- if (ctl$history && group) matrix(NA_real_, ctl$max_iter, 7L) else NULL
    h <- 0L
    Ct <- Wk <- Wl <- vector("list", length(E))
    for (it in seq_len(ctl$max_iter)) {
      check <- it == 1L || it == ctl$max_iter || it %% ctl$check_every == 0L
      shift <- if (group) 2 * mu else mu
      if (!identical(shift, last_shift)) {
        den <- lapply(z$D, function(D) D + shift); last_shift <- shift
      }
      rp2 <- rd2 <- nc2 <- nz2 <- ny2 <- 0
      for (e in E) {
        T <- if (group) (s$Gk[[e]] - s$Vk[[e]] + s$Gl[[e]] - s$Vl[[e]]) / 2 else
          s$Z[[e]] - s$H[[e]]
        Pt <- egcar_project(z, e, T)
        Ct[[e]] <- (z$St[[e]] + shift * Pt) / den[[e]]
        s$C[[e]] <- if (z$full[[e]]) egcar_lift(z, e, Ct[[e]]) else
          T + z$remainder[[e]] / shift + egcar_lift(z, e, Ct[[e]] - Pt)
        if (!group) {
          W <- s$C[[e]] + s$H[[e]]
          Zn <- soft_threshold(W, ctl$penalty / mu)
          Hn <- W - Zn
          if (check) {
            rp2 <- rp2 + sum((s$C[[e]] - Zn)^2)
            rd2 <- rd2 + sum((Zn - s$Z[[e]])^2)
            nc2 <- nc2 + sum(s$C[[e]]^2); nz2 <- nz2 + sum(Zn^2)
            ny2 <- ny2 + sum(Hn^2)
          }
          s$Z[[e]] <- Zn; s$H[[e]] <- Hn
        } else {
          Wk[[e]] <- s$C[[e]] + s$Vk[[e]]
          Wl[[e]] <- s$C[[e]] + s$Vl[[e]]
        }
      }
      if (group) {
        norms <- egcar_group_norms(z, Wk, Wl)
        threshold <- ctl$penalty / mu
        mult <- if (threshold <= 0) lapply(z$p_list, function(pk) rep.int(1, pk)) else
          lapply(norms, function(nr) pmax(0, 1 - threshold / pmax(nr, .Machine$double.eps)))
        for (e in E) {
          k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
          Gkn <- Wk[[e]] * mult[[k]]
          # sweep scales columns; no transpose of the large edge is required.
          Gln <- sweep(Wl[[e]], 2L, mult[[l]], "*")
          Vkn <- Wk[[e]] - Gkn; Vln <- Wl[[e]] - Gln
          if (check) {
            rp2 <- rp2 + sum((s$C[[e]] - Gkn)^2) + sum((s$C[[e]] - Gln)^2)
            rd2 <- rd2 + sum((Gkn - s$Gk[[e]] + Gln - s$Gl[[e]])^2)
            nc2 <- nc2 + 2 * sum(s$C[[e]]^2)
            nz2 <- nz2 + sum(Gkn^2) + sum(Gln^2)
            ny2 <- ny2 + sum((Vkn + Vln)^2)
          }
          s$Gk[[e]] <- Gkn; s$Gl[[e]] <- Gln
          s$Vk[[e]] <- Vkn; s$Vl[[e]] <- Vln
        }
      }
      if (check) {
        rp <- sqrt(rp2); rd <- mu * sqrt(rd2)
        ep <- sqrt(if (group) 2 * z$q else z$q) * ctl$abs_tol +
          ctl$rel_tol * max(sqrt(nc2), sqrt(nz2))
        ed <- sqrt(z$q) * ctl$abs_tol + ctl$rel_tol * mu * sqrt(ny2)
        if (ctl$history && group) {
          objective <- 0
          for (e in E) objective <- objective +
            0.5 * sum(z$D[[e]] * Ct[[e]]^2) - sum(z$S[[e]] * s$C[[e]])
          objective <- objective + ctl$penalty * sum(unlist(egcar_group_norms(z, s$C), use.names = FALSE))
          h <- h + 1L; history[h, ] <- c(it, objective, rp, rd, ep, ed, mu)
        }
        if (!is.finite(rp) || !is.finite(rd)) { warning("ADMM produced a non-finite residual."); break }
        if (rp <= ep && rd <= ed) { converged <- TRUE; break }
        if (ctl$adaptive && it %% ctl$adapt_every == 0L) {
          factor <- if (rp > ctl$balance_ratio * max(rd, .Machine$double.eps)) ctl$scale_factor else
            if (rd > ctl$balance_ratio * max(rp, .Machine$double.eps)) 1 / ctl$scale_factor else 1
          if (factor != 1) {
            mu <- mu * factor
            if (group) {
              s$Vk <- lapply(s$Vk, function(A) A / factor)
              s$Vl <- lapply(s$Vl, function(A) A / factor)
            } else s$H <- lapply(s$H, function(A) A / factor)
          }
        }
        if (verbose && (it == 1L || it %% 100L == 0L)) cat(sprintf(
          "  iter=%d primal=%.3e (%.3e) dual=%.3e (%.3e)\n", it, rp, ep, rd, ed))
      }
    }
    list(state = s, converged = converged, iterations = it, primal = rp, dual = rd,
         eps_primal = ep, eps_dual = ed, mu = mu,
         history = if (is.null(history)) matrix(numeric(), 0L, 7L) else history[seq_len(h), , drop = FALSE])
  }

  egcar_controls <- function(penalty, mu, max_iter, abs_tol, rel_tol, adaptive_mu,
                             balance_ratio, scale_factor, adapt_every, check_every, keep_history) {
    scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!scalar(penalty) || penalty < 0 || !scalar(mu) || mu <= 0)
      stop("EGCAR requires a finite nonnegative penalty and a finite positive ADMM parameter.")
    if (!scalar(abs_tol) || abs_tol < 0 || !scalar(rel_tol) || rel_tol < 0)
      stop("EGCAR tolerances must be finite and nonnegative.")
    for (x in list(max_iter, adapt_every, check_every))
      if (!scalar(x) || x < 1 || x != floor(x)) stop("EGCAR iteration counts must be positive integers.")
    if (!scalar(balance_ratio) || balance_ratio <= 0 || !scalar(scale_factor) || scale_factor <= 1)
      stop("Invalid residual-balancing parameters.")
    list(penalty = penalty, mu = mu, max_iter = as.integer(max_iter), abs_tol = abs_tol,
      rel_tol = rel_tol, adaptive = isTRUE(adaptive_mu), balance_ratio = balance_ratio,
      scale_factor = scale_factor, adapt_every = as.integer(adapt_every),
      check_every = as.integer(check_every), history = isTRUE(keep_history))
  }

  .egcar_check_solver_state <- function(state, context, group, backend) {
    fields <- if (group) c("C", "Gk", "Gl", "Vk", "Vl") else c("C", "Z", "H")
    hint <- if (identical(backend, "cpp")) paste0(
      " Reinstall the patched egcar source package and restart R, including CV workers.") else ""
    fail <- function(message) stop("EGCAR ", backend, " matrix interface: ", message,
                                    hint, call. = FALSE)
    if (!is.list(state) || anyDuplicated(names(state)) ||
        !all(fields %in% names(state)))
      fail("the solver did not return the required named state lists.")
    count <- length(context$edge_k)
    for (nm in fields) {
      blocks <- state[[nm]]
      if (!is.list(blocks) || length(blocks) != count)
        fail(paste0("state$", nm, " must contain ", count, " edge matrices."))
      for (j in seq_len(count)) {
        expected <- as.integer(c(context$p_list[[context$edge_k[[j]]]],
                                 context$p_list[[context$edge_l[[j]]]]))
        A <- blocks[[j]]
        if (!is.matrix(A) || !is.numeric(A) || !identical(dim(A), expected)) {
          observed <- if (is.null(dim(A))) paste0(
            "a dimensionless ", typeof(A), " object of length ", length(A)) else paste0(
            "a ", typeof(A), " array with dimensions ", paste(dim(A), collapse = " x "))
          fail(paste0("state$", nm, "[[", j, "]] (edge ",
            context$edge_k[[j]], "_", context$edge_l[[j]], ") must be a numeric ",
            paste(expected, collapse = " x "), " matrix; got ", observed, "."))
        }
      }
    }
    invisible(TRUE)
  }

  egcar_run_solver <- function(prep, ctl, init, group, verbose = FALSE) {
    z <- egcar_get_context(prep)
    s <- egcar_initial_state(prep, z, init, group)
    fun <- if (EGCAR_BACKEND == "cpp" ||
      (EGCAR_BACKEND == "auto" && identical(EGCAR_MASTER_BACKEND, "cpp")))
      egcar_native_function(required = EGCAR_BACKEND == "cpp") else NULL
    raw <- if (is.function(fun)) fun(z, s, ctl, group, verbose) else
      egcar_solve_R(z, s, ctl, group, verbose)
    if (is.function(fun) && !identical(raw$matrix_api, 2L))
      stop("EGCAR native matrix API mismatch. Restart R and use the updated standalone script.",
           call. = FALSE)
    .egcar_check_solver_state(raw$state, z, group,
                              if (is.function(fun)) "cpp" else "r")
    raw$context <- z
    raw
  }

  fit_l11_admm <- function(prep, rho_e, mu = 1, max_iter = 2000L,
      abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
      scale_factor = 2, adapt_every = 10L, entry_zero_tol = 1e-10,
      init = NULL, verbose = FALSE, check_every = 1L) {
    if (EGCAR_BACKEND == "reference") return(fit_l11_admm_reference(
      prep, rho_e, mu, max_iter, abs_tol, rel_tol, adaptive_mu, balance_ratio,
      scale_factor, adapt_every, entry_zero_tol, init, verbose, check_every))
    ctl <- egcar_controls(rho_e, mu, max_iter, abs_tol, rel_tol, adaptive_mu,
      balance_ratio, scale_factor, adapt_every, check_every, FALSE)
    a <- egcar_run_solver(prep, ctl, init, FALSE, verbose)
    s <- a$state; keys <- a$context$keys
    for (nm in c("C", "Z", "H")) names(s[[nm]]) <- keys
    C_hat <- lapply(s$Z, function(A) { A[abs(A) < entry_zero_tol] <- 0; A })
    list(C_hat = C_hat, C = s$C, Z = s$Z, H = s$H,
      converged = a$converged, iterations = a$iterations, primal_residual = a$primal,
      dual_residual = a$dual, eps_primal = a$eps_primal, eps_dual = a$eps_dual,
      rho_e = rho_e, lambda_g = 0, mu_z = a$mu)
  }

  fit_l21_admm <- function(prep, lambda_g, mu_g = 1, max_iter = 2000L,
      abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
      scale_factor = 2, adapt_every = 10L, group_zero_tol = 1e-8,
      entry_zero_tol = 1e-10, init = NULL, keep_history = FALSE,
      verbose = FALSE, check_every = 1L) {
    if (EGCAR_BACKEND == "reference") return(fit_l21_admm_reference(
      prep, lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu, balance_ratio,
      scale_factor, adapt_every, group_zero_tol, entry_zero_tol, init,
      keep_history, verbose, check_every))
    ctl <- egcar_controls(lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu,
      balance_ratio, scale_factor, adapt_every, check_every, keep_history)
    a <- egcar_run_solver(prep, ctl, init, TRUE, verbose)
    s <- a$state; z <- a$context; names(s$C) <- z$keys
    G <- egcar_view_copies(z, s$Gk, s$Gl)
    V <- egcar_view_copies(z, s$Vk, s$Vl)
    active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
    C_hat <- lapply(s$C, function(A) { A[abs(A) < entry_zero_tol] <- 0; A })
    for (e in seq_along(z$edge_k)) {
      C_hat[[e]][!active_rows[[z$edge_k[[e]]]], ] <- 0
      C_hat[[e]][, !active_rows[[z$edge_l[[e]]]]] <- 0
    }
    history <- if (keep_history) {
      h <- as.data.frame(a$history)
      names(h) <- c("iter", "objective", "primal", "dual", "eps_primal", "eps_dual", "mu_g")
      h$iter <- as.integer(h$iter); h
    } else NULL
    list(C_hat = C_hat, C = s$C, G = G, V = V, active_rows = active_rows,
      converged = a$converged, iterations = a$iterations, primal_residual = a$primal,
      dual_residual = a$dual, eps_primal = a$eps_primal, eps_dual = a$eps_dual,
      rho_e = 0, lambda_g = lambda_g, mu_g = a$mu, history = history)
  }

  # Same localized symmetric normalization as loading_from_operator(), but exploit
  # the block-diagonal metric, selected rows, and (on larger supports) a partial
  # largest-ALGEBRAIC eigensolve. Never replace this by the two-view SVD of ccar3.
  egcar_loading_factors <- function(prep, selected, covariance_ridge) {
    key <- paste0("egcar-block:", sprintf("%.17g", covariance_ridge), ":", paste(selected, collapse = ","))
    cache <- prep$loading_factor_cache
    if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))
    idx <- prep$indices %||% make_block_indices(prep$p_list)
    local <- lapply(idx, function(ii) which(ii %in% selected))
    position <- lapply(seq_along(idx), function(k) match(idx[[k]][local[[k]]], selected))
    diag_selected <- unlist(lapply(seq_along(idx), function(k) diag(prep$S_kk[[k]])[local[[k]]]), use.names = FALSE)
    scale_diag <- mean(diag_selected)
    if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
    blocks <- lapply(seq_along(idx), function(k) {
      sk <- local[[k]]
      if (!length(sk)) return(list(half = matrix(0, 0L, 0L), inv_half = matrix(0, 0L, 0L)))
      ev <- eigen(symmetrize(prep$S_kk[[k]][sk, sk, drop = FALSE]), symmetric = TRUE)
      d <- pmax(ev$values + covariance_ridge * scale_diag, 1e-10)
      list(half = tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors),
           inv_half = tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors))
    })
    ans <- list(local = local, position = position, blocks = blocks)
    if (is.environment(cache) && length(cache) < LOADING_FACTOR_CACHE_MAX) assign(key, ans, envir = cache)
    ans
  }

  egcar_top_eigen <- function(A, rank, positive_tol = 1e-10) {
    n <- nrow(A)
    # Full eigen is usually cheaper for the small default p=45 experiment.
    use_partial <- EGCAR_PARTIAL_EIGEN && n >= EGCAR_PARTIAL_EIGEN_MIN &&
      4L * rank < n && requireNamespace("RSpectra", quietly = TRUE)
    if (use_partial) {
      ee <- tryCatch(suppressWarnings(RSpectra::eigs_sym(A, k = rank, which = "LA",
        opts = list(tol = 1e-12, maxitr = 2000L,
                    ncv = min(n, max(4L * rank + 1L, 30L)),
                    initvec = sin(seq_len(n)) + cos(sqrt(2) * seq_len(n))))),
        error = function(e) NULL)
      if (!is.null(ee) && ee$nconv == rank && length(ee$values) == rank &&
          all(is.finite(ee$values)) && all(is.finite(ee$vectors))) {
        ord <- order(ee$values, decreasing = TRUE)
        ee$values <- ee$values[ord]; ee$vectors <- ee$vectors[, ord, drop = FALSE]
        residual <- A %*% ee$vectors - sweep(ee$vectors, 2L, ee$values, "*")
        error <- sqrt(colSums(residual * residual))
        scale <- max(1, max(rowSums(abs(A))))
        orth_error <- max(abs(crossprod(ee$vectors) - diag(rank)))
        # Resolve numerical rank/positivity near the original threshold using
        # the original dense solver, rather than letting iterative error decide.
        away_from_cutoff <- abs(ee$values[[rank]] - positive_tol) >
          max(1e-8 * scale, 4 * max(error))
        if (max(error) <= 1e-9 * scale && orth_error <= 1e-8 && away_from_cutoff)
          return(list(values = ee$values, vectors = ee$vectors))
      }
    }
    ee <- eigen(A, symmetric = TRUE)
    list(values = ee$values[seq_len(rank)], vectors = ee$vectors[, seq_len(rank), drop = FALSE])
  }

  egcar_loading_from_operator <- function(prep, C, rank, row_threshold = 1e-4,
      covariance_ridge = 1e-4, require_positive = TRUE, positive_tol = 1e-10,
      keep_full_C = TRUE) {
    if (EGCAR_BACKEND == "reference") return(loading_from_operator(
      prep, C, rank, row_threshold, covariance_ridge, require_positive, positive_tol))
    z <- egcar_get_context(prep)
    norms <- egcar_group_norms(z, C)
    selected <- which(unlist(norms, use.names = FALSE) > row_threshold)
    if (length(selected) < rank) return(list(valid = FALSE, reason = "fewer selected rows than rank"))
    f <- egcar_loading_factors(prep, selected, covariance_ridge)
    Rsel <- matrix(0, length(selected), length(selected))
    for (e in seq_along(z$edge_k)) {
      k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
      sk <- f$local[[k]]; sl <- f$local[[l]]
      if (!length(sk) || !length(sl)) next
      B <- f$blocks[[k]]$half %*% C[[e]][sk, sl, drop = FALSE] %*% f$blocks[[l]]$half
      Rsel[f$position[[k]], f$position[[l]]] <- B
      Rsel[f$position[[l]], f$position[[k]]] <- t(B)
    }
    ee <- egcar_top_eigen(symmetrize(Rsel), rank, positive_tol)
    if (require_positive && ee$values[[rank]] <= positive_tol)
      return(list(valid = FALSE, reason = "fewer than rank positive eigenvalues"))
    U <- ee$vectors
    L <- matrix(0, prep$p, rank)
    for (k in seq_along(z$p_list)) {
      ii <- f$position[[k]]
      if (length(ii)) L[selected[ii], ] <- f$blocks[[k]]$inv_half %*% U[ii, , drop = FALSE]
    }
    list(valid = TRUE, L = L, U = U, selected = selected, eigenvalues = ee$values,
      generalized_eigenvalues = 1 + ee$values,
      C_full = if (keep_full_C) assemble_full_C(C, prep$p_list) else NULL)
  }

  # Build/load once before benchmark timers, then LOAD the same build per worker.
  # Worker-local option storage avoids exporting Rcpp external pointers via future.
  egcar_startup <- proc.time()[[3L]]
  EGCAR_MASTER_BACKEND <- if (EGCAR_BACKEND %in% c("auto", "cpp")) {
    dir.create(EGCAR_CPP_CACHE, recursive = TRUE, showWarnings = FALSE)
    existing <- if (file.exists(EGCAR_CPP_FILE)) paste(readLines(EGCAR_CPP_FILE, warn = FALSE), collapse = "\n") else ""
    if (!identical(existing, EGCAR_CPP_SOURCE)) writeLines(EGCAR_CPP_SOURCE, EGCAR_CPP_FILE)
    if (is.function(egcar_native_function(required = EGCAR_BACKEND == "cpp"))) "cpp" else "r"
  } else EGCAR_BACKEND
  EGCAR_WORKER_BACKENDS <- if (PARALLEL_CV && EGCAR_MASTER_BACKEND == "cpp") {
    unlist(parallel_map_candidates(seq_len(CV_WORKERS), function(i) {
      set_blas_threads_one()
      if (is.function(egcar_native_function(required = EGCAR_BACKEND == "cpp"))) "cpp" else "r"
    }), use.names = FALSE)
  } else rep(EGCAR_MASTER_BACKEND, CV_WORKERS)
  EGCAR_ACCELERATION_INFO <- list(
    requested_backend = EGCAR_BACKEND, master_backend = EGCAR_MASTER_BACKEND,
    worker_backends = EGCAR_WORKER_BACKENDS,
    startup_seconds = proc.time()[[3L]] - egcar_startup,
    spectral_reduction = "thin SVD for wide training views; cached eigenbases otherwise; no user truncation; full null-space correction",
    reference = "cran/ccar3 R/ecca.r; inspected 2026-09-10; version 0.1.2",
    partial_loading_eigen = EGCAR_PARTIAL_EIGEN,
    native_error = if (EGCAR_MASTER_BACKEND == "r" && EGCAR_BACKEND == "auto")
      getOption("egcar.native.runtime.v2.matrix.api")$error else NULL)
  saveRDS(EGCAR_ACCELERATION_INFO, file.path(OUT_DIR, "egcar_acceleration_info.rds"))
  cat(sprintf("EGCAR backend: %s; startup (not fit time): %.2fs; CV workers: %s\n",
    EGCAR_MASTER_BACKEND, EGCAR_ACCELERATION_INFO$startup_seconds,
    paste(EGCAR_WORKER_BACKENDS, collapse = ",")))
  if (EGCAR_BACKEND == "auto" && EGCAR_MASTER_BACKEND == "r")
    message("Using optimized R ADMM. For compiled acceleration install Rcpp and RcppArmadillo, with a C++ compiler.")



    if (definitions_only) return(environment())

  config <- list(
    p_list = P_LIST,
    n_grid = N_GRID,
    rank_grid = RANK_GRID,
    active_per_view = ACTIVE_PER_VIEW,
    toeplitz_rho = TOEPLITZ_RHO,
    signal = SIGNAL,
    rho_e_cv_grid = RHO_E_CV_GRID,
    lambda_g_cv_grid = LAMBDA_G_CV_GRID,
    rate_c_e = RATE_C_E,
    rate_c_g = RATE_C_G,
    n_folds = N_FOLDS,
    n_reps = N_REP,
    workers = N_WORKERS,
    workers_requested = N_WORKERS_REQUESTED,
    cv_workers = CV_WORKERS,
    cv_workers_by_method = CV_WORKERS_BY_METHOD,
    cv_worker_allocation = CV_WORKER_ALLOCATION,
    cv_parallelism = "all CV methods sequential; same shared fold-level worker pool",
    plot_method_colors = METHOD_COLORS,
    plot_line_type = "solid",
    plot_point_shape = 16L,
    row_threshold = ROW_THRESHOLD,
    covariance_ridge = COVARIANCE_RIDGE,
    max_iter_cv = MAX_ITER_CV,
    max_iter_final = MAX_ITER_FINAL,
    run_external_benchmarks = RUN_EXTERNAL_BENCHMARKS,
    external_benchmark_methods = EXTERNAL_BENCHMARK_METHODS,
    common_cv_loss = "negative of validation_score(L, shared_fold$validation), with identical ridge/eigenvalue floor",
    external_preprocessing = "training-mean centering; no variable or block scaling",
    external_parallelism = "methods sequential; parallel over shared folds",
    sgca_k_grid = SGCA_K_GRID,
    sgca_rho_grid = SGCA_RHO_GRID,
    sgca_rho_rule = "direct coefficients, not rate multipliers",
    sgca_lambda_grid = SGCA_LAMBDA_GRID,
    sgca_tgd_algorithm = "Gao--Ma Algorithm 1; scaled iterate W=sqrt(lambda)V; backtracking; no per-iteration normalization",
    sgca_eta = SGCA_ETA,
    sgca_ridge_b = SGCA_RIDGE_B,
    sgca_init_tol = SGCA_INIT_TOL,
    sgca_max_iter_init = SGCA_MAX_ITER_INIT,
    sgca_tgd_tol = SGCA_TGD_TOL,
    sgca_max_iter_tgd = SGCA_MAX_ITER_TGD,
    rgcca_tau_grid = RGCCA_TAU_GRID,
    rgcca_scheme = RGCCA_SCHEME,
    rgcca_tol = RGCCA_TOL,
    rgcca_max_iter = RGCCA_MAX_ITER,
    sgcca_sparsity_grid = SGCCA_SPARSITY_GRID,
    multicca_l1_grid = MULTICCA_L1_GRID,
    multicca_niter = MULTICCA_NITER,
    external_training_sign_alignment = ALIGN_EXTERNAL_BLOCK_SIGNS,
    fast_sgca_initializer = FAST_SGCA_INITIALIZER,
    multicca_matrix_backend = MULTICCA_BACKEND,
    loading_factor_cache_max = LOADING_FACTOR_CACHE_MAX,
    make_loading_plots = MAKE_LOADING_PLOTS,
    loading_plot_n = LOADING_PLOT_N,
    loading_plot_ranks = LOADING_PLOT_RANKS,
    loading_plot_reps = LOADING_PLOT_REPS,
    loading_methods_per_page = LOADING_METHODS_PER_PAGE,
    loading_visualization_normalization = "global Euclidean orthonormal basis; global Procrustes to truth for display only",
    smoke_test = SMOKE_TEST,
    package_versions = setNames(vapply(c("RGCCA", "PMA", "future", "future.apply", "ggplot2"),
      function(pkg) if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_,
      character(1L)), c("RGCCA", "PMA", "future", "future.apply", "ggplot2")),
    oracle1_max_iter = ORACLE1_MAX_ITER,
    oracle1_algorithm = "zero-penalty consensus ADMM on population covariance blocks",
    oracle1_penalty_rule = "rho_e=0 and lambda_g=0; no cross-validation"
  )

  config$egcar_acceleration <- EGCAR_ACCELERATION_INFO

  results <- data.frame()
  cv_results <- data.frame()
  egcar_cv_fold_results <- data.frame()
  external_cv_fold_results <- data.frame()
  fits_store <- list()
  populations <- list()

  for (rep_id in seq_len(N_REP)) {
    for (rank in RANK_GRID) {
      population_seed <- MASTER_SEED + 100000L * rep_id + 1000L * rank
      data_seed <- population_seed + 77L
      population <- make_population(
        p_list = P_LIST,
        rank = rank,
        active_per_view = ACTIVE_PER_VIEW,
        toeplitz_rho = TOEPLITZ_RHO,
        signal = SIGNAL,
        seed = population_seed
      )
      pop_tag <- paste0("rep", rep_id, "_r", rank)
      populations[[pop_tag]] <- population

      # Oracle1 uses the exact population covariance blocks and zero values for
      # both penalty coefficients.  It therefore does not depend on n and is
      # fitted only once for the current replication and rank.
      population_prep <- prepare_population_problem(population)
      oracle1 <- fit_oracle_population(
        population = population,
        rank = rank,
        prep = population_prep,
        keep_history = TRUE
      )
      if (!identical(oracle1$status, "ok")) {
        warning(sprintf(
          "Oracle1 zero-penalty population ADMM status for rep=%d, r=%d: %s%s",
          rep_id, rank, oracle1$status,
          if (is.na(oracle1$error)) "" else paste0("; ", oracle1$error)
        ))
      }

      # Nested samples: each larger n extends the same simulated realization.
      all_views <- simulate_views(population, max(N_GRID), data_seed)

      for (n in N_GRID) {
        cat(sprintf("\nrep=%d  r=%d  n=%d\n", rep_id, rank, n))
        views <- lapply(all_views, function(X) X[seq_len(n), , drop = FALSE])
        centered_full <- center_views(views)$views
        full_prep <- prepare_problem(centered_full)
        fold_seed <- MASTER_SEED + 1000000L * rep_id + 10000L * rank + n
        fold_id <- make_folds(n, N_FOLDS, fold_seed)
        fold_objects <- make_fold_objects(views, fold_id)
        p <- full_prep$p
        d_max <- max(p - full_prep$p_list)
        base_e <- sqrt(log(p) / n)
        base_g <- sqrt((d_max + log(p)) / n)

        point_rows <- list()

        # Oracle1: population covariance inputs and exactly zero coefficients.
        # The same n-independent fit is evaluated at every n for plotting.
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "Oracle1-population",
          C = oracle1$C,
          loading = oracle1$loading,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = oracle1$time, tune_time = 0,
          rho_e = 0, lambda_g = 0,
          converged = oracle1$converged,
          iterations = oracle1$iterations,
          status = oracle1$status,
          error_message = oracle1$error
        )

        # Oracle2: empirical fit restricted to the true supports.
        start <- proc.time()[[3L]]
        C_oracle2 <- fit_oracle_support(full_prep, population$active_local)
        oracle2_time <- proc.time()[[3L]] - start
        L_oracle2 <- loading_from_operator(
          full_prep, C_oracle2, rank,
          row_threshold = ROW_THRESHOLD,
          covariance_ridge = COVARIANCE_RIDGE
        )
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "Oracle2-support",
          C = C_oracle2,
          loading = L_oracle2,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = oracle2_time, tune_time = 0,
          converged = TRUE, iterations = 1L
        )

        # Entrywise-only EGCAR with rate-scaled penalty.
        start <- proc.time()[[3L]]
        fit_e_rate <- fit_estimator(
          full_prep,
          rho_e = RATE_C_E * base_e,
          lambda_g = 0,
          max_iter = MAX_ITER_FINAL
        )
        time_e_rate <- proc.time()[[3L]] - start
        L_e_rate <- egcar_loading_from_operator(
          full_prep, fit_e_rate$C_hat, rank,
          ROW_THRESHOLD, COVARIANCE_RIDGE
        )
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "EGCAR-L11-rate",
          C = fit_e_rate$C_hat,
          loading = L_e_rate,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = time_e_rate,
          rho_e = RATE_C_E * base_e,
          lambda_g = 0,
          c_e = RATE_C_E, c_g = 0,
          converged = fit_e_rate$converged,
          iterations = fit_e_rate$iterations
        )

        # Entrywise-only EGCAR with direct one-dimensional CV over rho_e.
        cv_e <- cross_validate_penalties(
          full_views = views,
          full_prep = full_prep,
          fold_objects = fold_objects,
          rank = rank,
          rho_e_grid = RHO_E_CV_GRID,
          lambda_g_grid = 0,
          method = "l11",
          max_iter_cv = MAX_ITER_CV,
          max_iter_final = MAX_ITER_FINAL
        )
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "EGCAR-L11-CV",
          C = cv_e$fit$C_hat,
          loading = cv_e$loading,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = cv_e$fit_time,
          tune_time = cv_e$tuning_time,
          rho_e = cv_e$rho_e,
          lambda_g = 0,
          converged = cv_e$fit$converged,
          iterations = cv_e$fit$iterations,
          status = cv_e$status, error_message = cv_e$error
        )

        cv_e_table <- annotate_cv_table(
          cv_e$cv_table, "EGCAR-L11-CV", rep_id, rank, n,
          cv_e$best$candidate %||% NA_integer_
        )
        cv_results <- bind_rows_fill(cv_results, cv_e_table)

        # L21-only EGCAR with the retained group-rate scaling.
        start <- proc.time()[[3L]]
        fit_g_rate <- fit_estimator(
          full_prep,
          rho_e = 0,
          lambda_g = RATE_C_G * base_g,
          max_iter = MAX_ITER_FINAL,
          keep_history = TRUE,
          l21_only = TRUE
        )
        time_g_rate <- proc.time()[[3L]] - start
        L_g_rate <- egcar_loading_from_operator(
          full_prep, fit_g_rate$C_hat, rank,
          ROW_THRESHOLD, COVARIANCE_RIDGE
        )
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "EGCAR-L21-rate",
          C = fit_g_rate$C_hat,
          loading = L_g_rate,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = time_g_rate,
          rho_e = 0,
          lambda_g = RATE_C_G * base_g,
          c_e = 0, c_g = RATE_C_G,
          converged = fit_g_rate$converged,
          iterations = fit_g_rate$iterations
        )

        # L21-only EGCAR: positive, direct one-dimensional CV over lambda_g.
        cv_g <- cross_validate_penalties(
          full_views = views,
          full_prep = full_prep,
          fold_objects = fold_objects,
          rank = rank,
          rho_e_grid = 0,
          lambda_g_grid = LAMBDA_G_CV_GRID,
          method = "l21",
          max_iter_cv = MAX_ITER_CV,
          max_iter_final = MAX_ITER_FINAL
        )
        point_rows[[length(point_rows) + 1L]] <- evaluate_method(
          method = "EGCAR-L21-CV",
          C = cv_g$fit$C_hat,
          loading = cv_g$loading,
          population = population,
          rank = rank, n = n, rep_id = rep_id,
          fit_time = cv_g$fit_time,
          tune_time = cv_g$tuning_time,
          rho_e = cv_g$rho_e,
          lambda_g = cv_g$lambda_g,
          converged = cv_g$fit$converged,
          iterations = cv_g$fit$iterations,
          status = cv_g$status, error_message = cv_g$error
        )

        cv_g_table <- annotate_cv_table(
          cv_g$cv_table, "EGCAR-L21-CV", rep_id, rank, n,
          cv_g$best$candidate %||% NA_integer_
        )
        cv_results <- bind_rows_fill(cv_results, cv_g_table)

        if (SAVE_CV_FOLD_RESULTS) {
          for (label in c("EGCAR-L11-CV", "EGCAR-L21-CV")) {
            cv_one <- if (label == "EGCAR-L11-CV") cv_e else cv_g
            tab <- cv_one$cv_fold_table
            if (nrow(tab)) {
              tab$rep <- rep_id; tab$rank <- rank; tab$n <- n; tab$method <- label
              egcar_cv_fold_results <- bind_rows_fill(egcar_cv_fold_results, tab)
            }
          }
          utils::write.csv(egcar_cv_fold_results,
            file.path(OUT_DIR, "egcar_cv_fold_results.csv"), row.names = FALSE)
        }

        # All four external methods use the SAME shared fold objects/loss.
        # They estimate loading spaces, not C, so C/support metrics remain NA.
        external <- run_external_benchmarks(
          views, P_LIST, rank, fold_seed, fold_objects, full_prep
        )
        for (label in EXTERNAL_BENCHMARK_METHODS) {
          one_benchmark <- external[[label]]
          ext_loading <- if (!is.null(one_benchmark$L)) {
            list(valid = TRUE, L = one_benchmark$L)
          } else {
            NULL
          }
          point_rows[[length(point_rows) + 1L]] <- evaluate_method(
            method = label,
            C = NULL,
            loading = ext_loading,
            population = population,
            rank = rank, n = n, rep_id = rep_id,
            fit_time = one_benchmark$fit_time,
            tune_time = one_benchmark$tuning_time,
            converged = one_benchmark$converged,
            iterations = one_benchmark$iterations,
            status = one_benchmark$status,
            error_message = one_benchmark$error
          )
          ext_cv <- annotate_cv_table(
            one_benchmark$cv_table, label, rep_id, rank, n,
            one_benchmark$best$candidate %||% NA_integer_
          )
          cv_results <- bind_rows_fill(cv_results, ext_cv)
          ext_folds <- one_benchmark$cv_fold_table
          if (SAVE_CV_FOLD_RESULTS && !is.null(ext_folds) && nrow(ext_folds) > 0L) {
            ext_folds$rep <- rep_id
            ext_folds$rank <- rank
            ext_folds$n <- n
            ext_folds$method <- label
            external_cv_fold_results <- bind_rows_fill(external_cv_fold_results, ext_folds)
          }
          if (!identical(one_benchmark$status, "ok")) {
            warning(sprintf(
              "%s failed for rep=%d, r=%d, n=%d: %s",
              label, rep_id, rank, n, one_benchmark$error
            ))
          }
        }

        results <- rbind(results, do.call(rbind, point_rows))

        fit_tag <- paste0("rep", rep_id, "_r", rank, "_n", n)
        point_loadings <- list(
          "Oracle1-population" = oracle1$loading,
          "Oracle2-support" = L_oracle2,
          "EGCAR-L11-rate" = L_e_rate,
          "EGCAR-L11-CV" = cv_e$loading,
          "EGCAR-L21-rate" = L_g_rate,
          "EGCAR-L21-CV" = cv_g$loading
        )
        for (label in EXTERNAL_BENCHMARK_METHODS) {
          # Single-bracket assignment retains the name even for a NULL matrix.
          point_loadings[label] <- list(external[[label]]$L)
        }
        point_loadings <- point_loadings[METHOD_ORDER]
        if (SAVE_COMPACT_LOADINGS) {
          compact_dir <- file.path(OUT_DIR, "compact_loadings")
          dir.create(compact_dir, recursive = TRUE, showWarnings = FALSE)
          compact_mats <- lapply(point_loadings, function(x) tryCatch(plot_loading_matrix(x), error = function(e) NULL))
          saveRDS(list(loadings = compact_mats, truth = population$Lstar,
                       p_list = population$p_list, active_global = population$active_global,
                       rank = rank, n = n, rep = rep_id, signal = population$signal),
                  file.path(compact_dir, paste0(fit_tag, ".rds")), compress = "xz")
        }
        if (SAVE_FITS) {
          fits_store[[fit_tag]] <- list(
            Cstar = population$Cstar,
            Oracle1 = oracle1,
            Oracle2 = C_oracle2,
            EGCAR_L11_rate = fit_e_rate,
            EGCAR_L11_CV = cv_e,
            EGCAR_L21_rate = fit_g_rate,
            EGCAR_L21_CV = cv_g,
            external_benchmarks = external,
            loadings = point_loadings,
            fold_id = fold_id
          )
        }

        save_checkpoint(results, cv_results, fits_store, populations, config,
                        external_cv_fold_results)
        # Post-fit graphics are OUTSIDE every fit/tuning timer. They never alter
        # the reported loading matrices, CV tables, selected parameters or errors.
        tryCatch(
          make_loading_visualizations(point_loadings, population, rank, n, rep_id,
                                      do.call(rbind, point_rows)),
          error = function(e) {
            warning(sprintf("Loading plots failed for %s: %s", fit_tag, conditionMessage(e)))
            plot_dir <- file.path(OUT_DIR, "loading_visualizations", fit_tag)
            dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
            writeLines(conditionMessage(e), file.path(plot_dir, "plot_error.txt"))
          }
        )
        cat(sprintf("  selected L11 coefficient: rho_e=%g; status=%s\n",
          cv_e$rho_e, cv_e$status))
        cat(sprintf("  selected L21 coefficient: lambda_g=%g; status=%s\n",
          cv_g$lambda_g, cv_g$status))
        rm(views, centered_full, full_prep, fold_id, fold_objects,
           C_oracle2, L_oracle2, fit_e_rate, L_e_rate, cv_e,
           fit_g_rate, L_g_rate, cv_g, external, point_rows, point_loadings)
        gc(verbose = FALSE)
      }
    }
  }

  make_all_plots(results)
  utils::capture.output(utils::sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
  cat("\nCompleted. Main output file:\n")
  cat(file.path(OUT_DIR, "egcar_simulation.rds"), "\n")
  if (MAKE_LOADING_PLOTS) cat("Loading comparisons:", file.path(OUT_DIR, "loading_visualizations"), "\n")

  invisible(list(results = results, cv_results = cv_results,
    egcar_cv_fold_results = egcar_cv_fold_results,
    external_cv_fold_results = external_cv_fold_results,
    fits = fits_store, populations = populations, config = config,
    checkpoint = file.path(OUT_DIR, "egcar_simulation.rds")))

}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_local_egcar_experiments(
    output_dir = if (length(args) >= 1L) args[[1L]] else "egcar_local_outputs",
    workers = if (length(args) >= 2L) as.integer(args[[2L]]) else 1L,
    n_reps = if (length(args) >= 3L) as.integer(args[[3L]]) else 1L)
}

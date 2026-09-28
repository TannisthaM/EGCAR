# One wall-clock budget for SGCA preparation, the complete CV grid, and refit.
# Supervision is outside the numerical process so BLAS/LAPACK need not poll R.

.egcar_run_with_deadline <- function(fun, args, seconds) {
  .egcar_optional_limit(seconds, "seconds")
  if (!is.finite(seconds))
    return(list(timed_out = FALSE, result = do.call(fun, args), elapsed = NA_real_))
  start <- proc.time()[[3L]]
  logs <- c(tempfile("sgca-out-"), tempfile("sgca-err-"))
  child <- NULL
  on.exit({
    if (!is.null(child) && child$is_alive()) {
      try(child$kill_tree(), silent = TRUE)
      if (child$is_alive()) try(child$kill(), silent = TRUE)
    }
    unlink(logs)
  }, add = TRUE)
  child <- callr::r_bg(fun, args = args, libpath = .libPaths(),
    stdout = logs[[1L]], stderr = logs[[2L]], supervise = TRUE,
    system_profile = FALSE, user_profile = FALSE)
  repeat {
    elapsed <- proc.time()[[3L]] - start
    # Test the deadline before retrieving a late result. Completion after the
    # budget cannot be relabeled successful merely because a result exists.
    if (elapsed >= seconds) {
      try(child$kill_tree(), silent = TRUE)
      if (child$is_alive()) child$kill()
      return(list(timed_out = TRUE, result = NULL,
                  elapsed = proc.time()[[3L]] - start))
    }
    if (!child$is_alive()) {
      result <- child$get_result() # Solver errors stay errors, never timeouts.
      return(list(timed_out = FALSE, result = result,
                  elapsed = proc.time()[[3L]] - start))
    }
    child$wait(timeout = max(1L, min(100L, ceiling(1000 * (seconds - elapsed)))))
  }
}

sgca_budget_phase <- function(phase, tuning_time = NULL) {
  path <- get0("SGCA_BUDGET_PROGRESS", inherits = TRUE, ifnotfound = NULL)
  if (is.null(path)) return(invisible(NULL))
  # Only the coordinator writes this compact checkpoint, once at CV start
  # and once before refitting. Iterations do no file I/O.
  temp <- paste0(path, ".new")
  saveRDS(list(phase = phase, tuning_time = tuning_time), temp)
  if (!file.rename(temp, path)) {
    unlink(path)
    if (!file.rename(temp, path)) unlink(temp)
  }
  invisible(NULL)
}

sgca_common_cv <- function(
    full_views, full_prep, fold_objects, rank,
    k_grid = SGCA_K_GRID, rho_grid = SGCA_RHO_GRID,
    lambda_grid = SGCA_LAMBDA_GRID, seed = 1L,
    parallel_folds = PARALLEL_CV) {
  args <- list(full_views = full_views, full_prep = full_prep,
    fold_objects = fold_objects, rank = rank, k_grid = k_grid,
    rho_grid = rho_grid, lambda_grid = lambda_grid,
    seed = seed, parallel_folds = parallel_folds)
  if (!is.finite(SGCA_TIME_LIMIT)) {
    out <- do.call(sgca_common_cv_run, args)
    out$timed_out <- FALSE
    out$time_limit <- Inf
    out$time_limit_scope <- "SGCA CV and full-sample refit"
    return(out)
  }
  keys <- c("SGCA_ETA", "SGCA_RIDGE_B", "SGCA_INIT_TOL", "SGCA_MAX_ITER_INIT",
    "SGCA_TGD_TOL", "SGCA_MAX_ITER_TGD", "FAST_SGCA_INITIALIZER",
    "COVARIANCE_RIDGE", "ALIGN_EXTERNAL_BLOCK_SIGNS", "CV_WORKERS",
    "PARALLEL_CV", "BLAS_THREADS")
  settings <- mget(keys, envir = environment(), inherits = TRUE)
  settings$RETAIN_CV_FOLD_TABLES <- get0("RETAIN_CV_FOLD_TABLES", inherits = TRUE, ifnotfound = TRUE)
  settings$RETAIN_BENCHMARK_FITS <- get0("RETAIN_BENCHMARK_FITS", inherits = TRUE, ifnotfound = TRUE)
  settings$PARALLEL_CV <- isTRUE(parallel_folds) && isTRUE(settings$PARALLEL_CV)
  if (!settings$PARALLEL_CV) settings$CV_WORKERS <- 1L
  progress <- tempfile("sgca-phase-", fileext = ".rds")
  on.exit(unlink(c(progress, paste0(progress, ".new"))), add = TRUE)
  future_options <- options()[intersect(names(options()),
    c("future.globals.maxSize", "parallelly.fork.enable", "future.fork.enable"))]
  run <- .egcar_run_with_deadline(
    function(args, settings, progress, future_options, allow_fork) {
      ns <- asNamespace("egcar")
      e <- get(".egcar_engine", ns)()
      list2env(settings, envir = e)
      e$SGCA_BUDGET_PROGRESS <- progress
      options(future_options)
      if (!allow_fork) Sys.setenv(EGCAR_MULTICORE = "0")
      get(".egcar_with_threads", ns)(settings$BLAS_THREADS,
        get(".egcar_with_workers", ns)(settings$CV_WORKERS,
          do.call(e$sgca_common_cv_run, args)))
    }, list(args = args, settings = settings, progress = progress,
            future_options = future_options,
            allow_fork = .egcar_worker_plan_is_multicore()), SGCA_TIME_LIMIT)
  if (!run$timed_out) {
    out <- run$result
    # Include supervision/startup/transfer overhead in reported method time.
    out$supervision_overhead <- max(0, run$elapsed - out$fit_time - out$tuning_time)
    out$tuning_time <- out$tuning_time + out$supervision_overhead
    out$time <- out$fit_time + out$tuning_time
  } else {
    phase <- if (file.exists(progress)) tryCatch(readRDS(progress),
      error = function(e) list(phase = "startup_or_cv")) else list(phase = "startup_or_cv")
    tuning <- if (identical(phase$phase, "refit")) phase$tuning_time else run$elapsed
    tuning <- min(run$elapsed, tuning)
    out <- list(L = NULL, loading = NULL, fit_full = NULL, best = NULL,
      cv_table = data.frame(), cv_fold_table = data.frame(),
      fit_time = run$elapsed - tuning, tuning_time = tuning, time = run$elapsed,
      status = "time_limit", converged = FALSE, iterations = NA_integer_,
      timeout_phase = phase$phase,
      error = sprintf(paste0("SGCA did not converge in real time: the complete CV and refit ",
        "did not finish within the %.6g-hour wall-clock budget (%.6g seconds)."),
        SGCA_TIME_LIMIT / 3600, SGCA_TIME_LIMIT))
  }
  out$timed_out <- run$timed_out
  out$time_limit <- SGCA_TIME_LIMIT
  out$time_limit_scope <- "SGCA CV and full-sample refit"
  out$wall_time <- run$elapsed
  out
}

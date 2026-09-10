# Full local study using installed package functions.
run_egcar_experiments <- function(
    output_dir = "egcar_package_outputs", workers = 1L, n_reps = 1L,
    config = egcar_experiment_config(), backend = c("cpp", "R", "reference"),
    smoke_test = FALSE) {
  backend <- match.arg(backend)
  .egcar_scalar(workers, "workers", 1, integer = TRUE)
  .egcar_scalar(n_reps, "n_reps", 1, integer = TRUE)
  .egcar_flag(smoke_test, "smoke_test")
  config <- validate_egcar_experiment_config(config)
  if (smoke_test) { config <- egcar_smoke_config(config); n_reps <- 1L }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create output_dir.")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  requested_workers <- as.integer(workers)
  workers <- as.integer(min(workers, config$n_folds))
  if (workers > 1L && (!requireNamespace("future", quietly = TRUE) ||
                       !requireNamespace("future.apply", quietly = TRUE))) {
    message("Parallel packages missing: all CV methods use one worker.")
    workers <- 1L
  }
  ctl <- egcar_control(backend = backend, max_iter = config$max_iter_final,
    max_iter_cv = config$max_iter_cv, mu = config$mu_z, abs_tol = config$abs_tol,
    rel_tol = config$rel_tol, adaptive_mu = config$adaptive_mu,
    check_every_cv = config$check_every_admm, row_threshold = config$row_threshold,
    covariance_ridge = config$covariance_ridge, group_zero_tol = config$group_zero_tol,
    entry_zero_tol = config$entry_zero_tol, loading_cache_max = config$loading_factor_cache_max)
  bm <- benchmark_control(sgca_eta = config$sgca_eta, sgca_ridge = config$sgca_ridge_b,
    sgca_init_tol = config$sgca_init_tol, sgca_init_max_iter = config$sgca_max_iter_init,
    sgca_tgd_tol = config$sgca_tgd_tol, sgca_tgd_max_iter = config$sgca_max_iter_tgd,
    fast_sgca_initializer = config$fast_sgca_initializer, rgcca_scheme = config$rgcca_scheme,
    rgcca_tol = config$rgcca_tol, rgcca_max_iter = config$rgcca_max_iter,
    multicca_niter = config$multicca_niter, multicca_backend = config$multicca_backend,
    align_signs = config$align_external_block_signs)
  e <- .egcar_engine(ctl, bm, workers)
  for (nm in c("evaluate_method", "fit_estimator", "cross_validate_penalties", "annotate_cv_table", "run_external_benchmarks", "summarize_metric", "save_metric_plot", "make_plot_set", "make_all_plots", "loading_euclidean_basis", "plot_loading_matrix", "save_multipage_loading_pdf", "make_loading_visualizations", "save_checkpoint")) {
    fn <- get(nm, envir = environment(run_egcar_experiments), inherits = FALSE)
    environment(fn) <- e
    assign(nm, fn, envir = e)
  }
  list2env(setNames(unclass(config), toupper(names(config))), envir = e)
  e$N_REP <- as.integer(n_reps); e$SMOKE_TEST <- smoke_test
  e$OUT_DIR <- output_dir; e$N_WORKERS_REQUESTED <- requested_workers
  e$N_WORKERS <- e$CV_WORKERS <- workers; e$PARALLEL_CV <- workers > 1L
  e$EGCAR_BACKEND <- e$EGCAR_MASTER_BACKEND <- tolower(backend)
  e$EGCAR_PARTIAL_EIGEN <- !identical(Sys.getenv("EGCAR_PARTIAL_EIGEN", "1"), "0")
  # Each worker loads the installed DLL by namespace; no external pointer is exported.
  .egcar_with_seed(NULL, .egcar_with_threads(1L, .egcar_with_workers(workers, {
    old_options <- options(stringsAsFactors = FALSE)
    on.exit(options(old_options), add = TRUE)
    old_max <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_max), add = TRUE)
    options(future.globals.maxSize = max(4 * 1024^3, old_max %||% 0))
    start <- proc.time()[[3L]]
    worker_backends <- if (workers > 1L) unlist(e$parallel_map_candidates(seq_len(workers),
      function(i) { loadNamespace("egcar"); tolower(backend) }), use.names = FALSE) else tolower(backend)
    e$EGCAR_ACCELERATION_INFO <- list(requested_backend = tolower(backend),
      master_backend = tolower(backend), worker_backends = worker_backends,
      startup_seconds = proc.time()[[3L]] - start,
      spectral_reduction = "thin SVD/cached eigenbases, full null-space correction",
      reference = "cran/ccar3 R/ecca.r computational techniques",
      partial_loading_eigen = e$EGCAR_PARTIAL_EIGEN,
      package_version = as.character(utils::packageVersion("egcar")))
    saveRDS(e$EGCAR_ACCELERATION_INFO, file.path(output_dir, "egcar_acceleration_info.rds"))
    evalq({
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
    }, envir = e)
    fun <- .egcar_experiment_loop
    environment(fun) <- e
    fun()
  })))
}

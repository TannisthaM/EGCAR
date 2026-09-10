# Canonical loop shared with the standalone runner.
.egcar_experiment_loop <- function() {
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

        for (label in c("EGCAR-L11-CV", "EGCAR-L21-CV")) {
          cv_one <- if (label == "EGCAR-L11-CV") cv_e else cv_g
          tab <- cv_one$cv_fold_table
          tab$rep <- rep_id; tab$rank <- rank; tab$n <- n; tab$method <- label
          egcar_cv_fold_results <- bind_rows_fill(egcar_cv_fold_results, tab)
        }
        utils::write.csv(egcar_cv_fold_results,
          file.path(OUT_DIR, "egcar_cv_fold_results.csv"), row.names = FALSE)

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
          if (!is.null(ext_folds) && nrow(ext_folds) > 0L) {
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

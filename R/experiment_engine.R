# Local-experiment adapters and plotting helpers; unchanged comparison CV routines live in alt_*.R.

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

fit_estimator <- function(prep, rho_e, lambda_g, max_iter, keep_history = FALSE,
                           init = NULL, check_every = 1L, l21_only = FALSE,
                           keep_state = TRUE) {
  if (isTRUE(l21_only)) {
    if (length(rho_e) != 1L || !is.finite(rho_e) || rho_e != 0)
      stop("L21-only EGCAR requires rho_e = 0.")
    return(fit_l21_admm(prep, lambda_g, mu_g = MU_G, max_iter = max_iter,
      abs_tol = ABS_TOL, rel_tol = REL_TOL, adaptive_mu = ADAPTIVE_MU,
      group_zero_tol = GROUP_ZERO_TOL, entry_zero_tol = ENTRY_ZERO_TOL,
      init = init, check_every = check_every, keep_history = keep_history,
      keep_state = keep_state))
  }
  if (length(lambda_g) != 1L || !is.finite(lambda_g) || lambda_g != 0)
    stop("L11-only EGCAR requires lambda_g = 0; combined penalties are not supported.")
  fit_l11_admm(prep, rho_e, mu = MU_Z, max_iter = max_iter,
    abs_tol = ABS_TOL, rel_tol = REL_TOL, adaptive_mu = ADAPTIVE_MU,
    entry_zero_tol = ENTRY_ZERO_TOL, init = init, check_every = check_every,
    keep_state = keep_state)
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
  per_fold <- parallel_map_candidates(fold_objects, function(fo) {
    set_blas_threads_one()
    f <- fo$fold
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
    return(list(fit = missing_fit, loading = NULL, cv_table = table,
      cv_fold_table = if (isTRUE(SAVE_CV_FOLD_RESULTS)) fold_table else data.frame(),
      best = NULL, rho_e = if (group) 0 else NA_real_, lambda_g = if (group) NA_real_ else 0,
      fit_time = 0, tuning_time = tuning_time, status = "no_valid_cv",
      error = "No positive CV candidate gave a finite requested-rank score on every fold."))
  }
  best <- table[best_index, , drop = FALSE]
  start <- proc.time()[[3L]]
  fit <- tryCatch(fit_estimator(full_prep, best$rho_e[[1L]], best$lambda_g[[1L]],
    max_iter_final, keep_history = FALSE, l21_only = group,
    keep_state = isTRUE(SAVE_FITS)),
    error = function(e) list(error = conditionMessage(e)))
  fit_time <- proc.time()[[3L]] - start
  if (is.null(fit$C_hat)) return(list(fit = missing_fit, loading = NULL,
    cv_table = table,
    cv_fold_table = if (isTRUE(SAVE_CV_FOLD_RESULTS)) fold_table else data.frame(),
    best = best, rho_e = best$rho_e[[1L]],
    lambda_g = best$lambda_g[[1L]], fit_time = fit_time, tuning_time = tuning_time,
    status = "refit_error", error = fit$error))
  loading <- tryCatch(egcar_loading_from_operator(full_prep, fit$C_hat, rank,
    ROW_THRESHOLD, COVARIANCE_RIDGE, require_positive = TRUE, keep_full_C = FALSE),
    error = function(e) list(valid = FALSE, reason = conditionMessage(e)))
  status <- if (!isTRUE(loading$valid)) "invalid_loading" else
    if (!isTRUE(fit$converged)) "not_converged" else "ok"
  list(fit = fit, loading = loading, cv_table = table,
    cv_fold_table = if (isTRUE(SAVE_CV_FOLD_RESULTS)) fold_table else data.frame(),
    best = best, rho_e = best$rho_e[[1L]], lambda_g = best$lambda_g[[1L]],
    tuning_time = tuning_time, fit_time = fit_time, status = status,
    error = if (!isTRUE(loading$valid)) loading$reason else
      if (!isTRUE(fit$converged)) "ADMM did not satisfy its stopping rule." else NA_character_)
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
    rotation = diag(rank), importance = rowSums(Qtrue * Qtrue),
    status = "population truth",
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
      projector_error <- sqrt(max(0, 2 * rank - 2 * sum(crossprod(Q, Qtrue)^2)))
      list(valid = TRUE, raw = L, Q = Q, aligned = aligned,
           rotation = rotation, importance = rowSums(Q * Q),
           status = status, converged = conv,
           alignment_error = frob(aligned - Qtrue),
           projector_error = projector_error, error = NA_character_)
    }, error = function(e) list(
      valid = FALSE, raw = tryCatch(plot_loading_matrix(loadings[[label]]),
                                    error = function(e) NULL),
      Q = matrix(NA_real_, p, rank), aligned = matrix(NA_real_, p, rank),
      rotation = NULL, importance = rep(NA_real_, p),
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
  if (isTRUE(SAVE_LOADING_DATA)) {
    saveRDS(list(raw_loadings = loadings, population = population, records = records,
                 rank = rank, n = n, rep_id = rep_id, fit_results = fit_results,
                 normalization = "global Euclidean orthonormal basis",
                 alignment = "one global orthogonal Procrustes rotation for display only"),
            file.path(output_dir, "loading_plot_data.rds"), compress = "xz")
  }
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
  if (isTRUE(SAVE_LOADING_DATA)) {
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
    if (isTRUE(SAVE_LOADING_DATA))
      "Raw, normalized, and aligned values are in loading_coefficients.csv and loading_plot_data.rds."
    else "Full coefficient/RDS plot data were not stored (SAVE_LOADING_DATA=FALSE).",
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

  difference_limit <- max(vapply(methods, function(label) {
    z <- records[[label]]
    if (!isTRUE(z$valid)) return(NA_real_)
    finite_max(tcrossprod(z$Q) - Ptrue, fallback = 0)
  }, numeric(1L)), na.rm = TRUE)
  if (!is.finite(difference_limit) || difference_limit <= 1e-12) difference_limit <- 1
  projection_groups <- split(methods, ceiling(seq_along(methods) / 4L))
  boundaries <- head(cumsum(pp), -1L) + 0.5
  view_centers <- cumsum(pp) - (pp - 1) / 2
  projection_pages <- lapply(projection_groups, function(group) {
    d <- do.call(rbind, lapply(group, function(label) {
      ij <- expand.grid(row_variable = seq_len(p), column_variable = seq_len(p), KEEP.OUT.ATTRS = FALSE)
      ij$method <- label
      z <- records[[label]]
      ij$difference <- if (isTRUE(z$valid)) as.vector(tcrossprod(z$Q) - Ptrue) else rep(NA_real_, p * p)
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

save_compact_loadings <- function(loadings, population, rank, n, rep_id, fit_results) {
  if (!isTRUE(SAVE_COMPACT_LOADINGS)) return(invisible(NULL))
  dir <- file.path(OUT_DIR, "compact_loadings")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  mats <- lapply(loadings, function(x) tryCatch(plot_loading_matrix(x), error = function(e) NULL))
  status <- if (is.data.frame(fit_results) && nrow(fit_results))
    fit_results[, intersect(c("method", "status", "converged"), names(fit_results)), drop = FALSE] else data.frame()
  payload <- list(loadings = mats, truth = population$Lstar,
    p_list = population$p_list, active_global = population$active_global,
    rank = rank, n = n, rep = rep_id, status = status)
  saveRDS(payload, file.path(dir, paste0("rep", rep_id, "_r", rank, "_n", n, ".rds")),
          compress = "xz")
  invisible(NULL)
}

save_checkpoint <- function(results, cv_results, fits, populations, config,
                            external_cv_fold_results = data.frame()) {
  checkpoint <- list(results = results, cv_results = cv_results, config = config,
    storage = list(full_fits = isTRUE(SAVE_FITS), cv_fold_results = isTRUE(SAVE_CV_FOLD_RESULTS),
                   loading_plot_data = isTRUE(SAVE_LOADING_DATA), compact_loadings = isTRUE(SAVE_COMPACT_LOADINGS)))
  if (isTRUE(SAVE_FITS)) {
    checkpoint$fits <- fits
    checkpoint$populations <- populations
  }
  if (isTRUE(SAVE_CV_FOLD_RESULTS)) checkpoint$external_cv_fold_results <- external_cv_fold_results
  tmp <- tempfile("egcar_checkpoint_", tmpdir = OUT_DIR, fileext = ".rds")
  saveRDS(checkpoint, tmp, compress = "xz")
  if (!file.rename(tmp, file.path(OUT_DIR, "egcar_simulation.rds"))) {
    if (!file.copy(tmp, file.path(OUT_DIR, "egcar_simulation.rds"), overwrite = TRUE))
      stop("Could not write checkpoint.")
    unlink(tmp)
  }
  utils::write.csv(results, file.path(OUT_DIR, "simulation_results.csv"), row.names = FALSE)
  if (nrow(cv_results) > 0L) {
    utils::write.csv(cv_results, file.path(OUT_DIR, "cv_grid_results.csv"), row.names = FALSE)
    if ("selected" %in% names(cv_results)) {
      chosen <- cv_results[which(cv_results$selected %in% TRUE), , drop = FALSE]
      utils::write.csv(chosen, file.path(OUT_DIR, "selected_cv_parameters.csv"), row.names = FALSE)
    }
  }
  if (isTRUE(SAVE_CV_FOLD_RESULTS) && nrow(external_cv_fold_results) > 0L) {
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

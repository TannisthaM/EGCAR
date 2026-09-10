#!/usr/bin/env Rscript
# Rscript verify_EGCAR_rewrites.R run_EGCAR_local.R [check_output] [workers]
# Run AFTER installing egcar 0.2.0 and the benchmark dependencies.
# Runs tiny experiments and compares standalone vs installed-package results.
verify_egcar_rewrites <- function(local_script, output_dir = "egcar_rewrite_checks", workers = 1L) {
  if (!requireNamespace("egcar", quietly = TRUE)) stop("Install egcar first.")
  if (utils::packageVersion("egcar") < "0.2.0") stop("Need egcar >= 0.2.0.")
  local_script <- normalizePath(local_script, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(local_script, envir = env)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  ns <- asNamespace("egcar")
  assert_close <- function(a, b, label, tol = 1e-6) {
    eq <- all.equal(a, b, tolerance = tol, check.attributes = FALSE)
    if (!isTRUE(eq)) stop(label, ": ", paste(eq, collapse = "; "))
  }
  assert_error <- function(expr) {
    z <- tryCatch(force(expr), error = identity)
    if (!inherits(z, "error")) stop("Expected a validation error.")
  }
  has_comparisons <- all(vapply(c("RGCCA", "PMA"), requireNamespace,
                                 logical(1L), quietly = TRUE))
  cfg <- egcar::egcar_experiment_config(p_list = c(5L, 5L, 5L),
    n_grid = 24L, rank_grid = 1L, active_per_view = 3L, n_folds = 2L,
    rho_e_cv_grid = c(0.001, 0.01), lambda_g_cv_grid = c(0.001, 0.01),
    max_iter_cv = 30L, max_iter_final = 60L, oracle1_max_iter = 60L,
    sgca_k_grid = c(6L, 12L), sgca_rho_grid = c(0.01, 0.1),
    sgca_lambda_grid = c(0.1, 1), sgca_max_iter_init = 40L, sgca_max_iter_tgd = 80L,
    rgcca_tau_grid = c(0.1, 1), sgcca_sparsity_grid = c(0.6, 1),
    multicca_l1_grid = c(1.5, 2), rgcca_max_iter = 80L, multicca_niter = 10L,
    run_external_benchmarks = has_comparisons,
    make_plots = requireNamespace("ggplot2", quietly = TRUE),
    make_loading_plots = requireNamespace("ggplot2", quietly = TRUE))
  assert_error(egcar::egcar_experiment_config(rho_e_cv_grid = c(0, 0.01)))
  assert_error(env$egcar_experiment_config(lambda_g_cv_grid = c(0, 0.01)))

  # Verify function bodies and formals, not environment addresses.
  defs <- env$run_local_egcar_experiments(file.path(output_dir, "definitions"),
    workers = 1L, config = cfg, backend = "R", definitions_only = TRUE)
  routine_names <- c("sgca_common_cv", "rgcca_family_common_cv", "rgcca_common_cv",
    "sgcca_common_cv", "multicca_common_cv", "cross_validate_loading_grid",
    "validation_score", "sgca_tgd_penalized", "multicca_gram_fit")
  for (nm in routine_names) {
    a <- get(nm, defs, inherits = FALSE); b <- get(nm, ns, inherits = FALSE)
    if (!identical(deparse(formals(a), width.cutoff = 500L), deparse(formals(b), width.cutoff = 500L)) ||
        !identical(deparse(body(a), width.cutoff = 500L), deparse(body(b), width.cutoff = 500L)))
      stop("Local/package function mismatch: ", nm)
  }
  native <- list()
  for (be in c("R", "cpp")) {
    a <- suppressWarnings(env$run_local_egcar_experiments(
      file.path(output_dir, paste0("local_", be)), workers = workers,
      n_reps = 1L, config = cfg, backend = be))
    b <- suppressWarnings(egcar::run_egcar_experiments(
      file.path(output_dir, paste0("package_", be)), workers = workers,
      n_reps = 1L, config = cfg, backend = be))
    compare_columns <- c("rep", "rank", "n", "method", "C_error", "C_relative_error",
      "subspace_euclidean", "subspace_sigma0", "support_precision", "support_recall",
      "support_fdp", "rho_e", "lambda_g", "converged", "iterations", "status")
    assert_close(a$results[, compare_columns], b$results[, compare_columns], paste(be, "results"))
    cv_cols <- intersect(c("candidate", "lambda", "rho_e", "lambda_g", "method", "mean_loss",
      "valid_folds", "converged_folds", "selected", "sgca_k", "sgca_rho", "sgca_lambda",
      "rgcca_tau", "sgcca_sparsity", "multicca_l1_bound"), names(a$cv_results))
    assert_close(a$cv_results[, cv_cols], b$cv_results[, cv_cols], paste(be, "CV tables"))
    for (tag in names(a$fits)) {
      if (!identical(a$fits[[tag]]$fold_id, b$fits[[tag]]$fold_id)) stop("Fold mismatch: ", tag)
      for (nm in c("EGCAR_L11_rate", "EGCAR_L21_rate")) {
        assert_close(a$fits[[tag]][[nm]]$C_hat, b$fits[[tag]][[nm]]$C_hat, paste(be, tag, nm))
      }
      for (nm in c("EGCAR_L11_CV", "EGCAR_L21_CV")) {
        assert_close(a$fits[[tag]][[nm]]$fit$C_hat,
                     b$fits[[tag]][[nm]]$fit$C_hat, paste(be, tag, nm))
      }
    }
    if (has_comparisons) {
      ext <- b$results[b$results$method %in% c("SGCA", "RGCCA", "SGCCA", "MultiCCA"), ]
      if (any(!is.finite(ext$subspace_euclidean)))
        stop("A comparison returned no valid loading; inspect benchmark_failures.csv.")
    }
    native[[be]] <- b
  }
  numeric_cols <- c("C_error", "C_relative_error", "subspace_euclidean", "subspace_sigma0")
  assert_close(native$R$results[, numeric_cols], native$cpp$results[, numeric_cols],
                "Accelerated R/native outputs", tol = 1e-5)
  saveRDS(list(passed = TRUE, workers = workers, comparisons_tested = has_comparisons,
              package_version = as.character(utils::packageVersion("egcar")),
              local_script = local_script), file.path(output_dir, "verification.rds"))
  utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  cat("Standalone/package checks passed for accelerated R and compiled backends.\n")
  if (!has_comparisons) cat("RGCCA/PMA missing: external comparisons were SKIPPED.\n")
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Supply the path to run_EGCAR_local.R.")
  verify_egcar_rewrites(args[[1L]],
    if (length(args) >= 2L) args[[2L]] else "egcar_rewrite_checks",
    if (length(args) >= 3L) as.integer(args[[3L]]) else 1L)
}

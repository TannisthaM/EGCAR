# Focused installed-package regression checks. Sourcing defines a function only.
# Run through examples/07_matrix_interface_check.R or call this function directly.
# No RGCCA, PMA, ggplot2 or testthat dependency is needed for these checks.
run_egcar_matrix_interface_check <- function(workers = 1L, check_cv = TRUE,
                                            verbose = TRUE) {
  if (!requireNamespace("egcar", quietly = TRUE)) stop("Install egcar first.")
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1 || workers != floor(workers)) stop("workers must be a positive integer.")
  if (!is.logical(check_cv) || length(check_cv) != 1L || is.na(check_cv))
    stop("check_cv must be TRUE or FALSE.")
  workers <- as.integer(min(workers, 5L))
  if (check_cv && workers > 1L) {
    missing <- c("future", "future.apply")[!vapply(c("future", "future.apply"),
      requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing)) stop("Install the parallel-check dependencies: ",
                              paste(missing, collapse = ", "))
  }
  engine <- getFromNamespace(".egcar_engine", "egcar")
  native <- getFromNamespace("egcar_native_solve", "egcar")
  state_check <- getFromNamespace(".egcar_check_solver_state", "egcar")
  scope_workers <- getFromNamespace(".egcar_with_workers", "egcar")
  ctl <- egcar::egcar_control(backend = "cpp", max_iter = 150L, max_iter_cv = 100L,
    abs_tol = 0, rel_tol = 0, adaptive_mu = TRUE, check_every = 5L,
    check_every_cv = 5L, partial_eigen = FALSE)
  e <- engine(ctl)
  say <- function(...) if (isTRUE(verbose)) message(...)
  assert <- function(ok, label) if (!isTRUE(ok)) stop(label, call. = FALSE)
  matrix_check <- function(A, nr, nc, label) {
    assert(is.matrix(A) && is.numeric(A) &&
      identical(dim(A), as.integer(c(nr, nc))) && all(is.finite(A)),
      paste0(label, " is not a finite ", nr, " x ", nc, " numeric matrix."))
  }
  agree <- function(a, b, label, tolerance = 1e-8) {
    assert(length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)),
           paste0(label, ": invalid comparison values."))
    error <- max(c(0, abs(as.numeric(a) - as.numeric(b))))
    assert(error <= tolerance * max(c(1, abs(as.numeric(b)))),
           paste0(label, ": numerical disagreement (maximum absolute error ", error, ")."))
    error
  }
  check_cv_result <- function(result, label) {
    notes <- unique(stats::na.omit(result$cv_fold_table$error_message))
    assert(!is.null(result$fit) && all(is.finite(result$cv_fold_table$score)),
      paste0(label, ": CV did not produce finite scores/refit. ",
             if (length(notes)) paste(head(notes, 5L), collapse = " | ") else
               result$error))
  }
  check_fit <- function(fit, prepared, rank, family, label, require_loading = TRUE) {
    z <- e$egcar_get_context(prepared$prep)
    for (j in seq_along(z$edge_k)) {
      nr <- z$p_list[[z$edge_k[[j]]]]; nc <- z$p_list[[z$edge_l[[j]]]]
      for (nm in if (family == "l11") c("C", "Z", "H") else "C")
        matrix_check(fit$solver[[nm]][[j]], nr, nc, paste(label, nm, j))
      matrix_check(fit$C[[j]], nr, nc, paste(label, "reported C", j))
    }
    if (family == "l21") {
      assert(is.null(fit$solver$Z) && is.null(fit$solver$H) && fit$rho_e == 0,
             paste0(label, ": L21 must not contain an entrywise penalty/state."))
      for (k in seq_along(z$p_list)) for (nm in c("G", "V"))
        matrix_check(fit$solver[[nm]][[k]], z$p_list[[k]],
          sum(z$p_list) - z$p_list[[k]], paste(label, nm, k))
    } else assert(fit$lambda_g == 0, paste0(label, ": L11 must not contain a group penalty."))
    if (!is.null(fit$L)) {
      matrix_check(fit$L, prepared$p, rank, paste(label, "stacked loading"))
      for (k in seq_along(prepared$p_list))
        matrix_check(fit$loadings[[k]], prepared$p_list[[k]], rank,
                     paste(label, "view loading", k))
    } else {
      assert(!require_loading && identical(fit$status, "invalid_loading"),
        paste0(label, ": no requested-rank loading; status=", fit$status,
               "; ", fit$error))
    }
    invisible(TRUE)
  }
  say("Checking egcar ", utils::packageVersion("egcar"), " from ", find.package("egcar"))

  # Direct .Call checks happen BEFORE public wrappers can reshape/name the state.
  cases <- list(
    list(label = "one_by_one", n = 24L, p = c(1L, 1L), rank = 1L, active = 1L),
    list(label = "row_column_singletons", n = 24L, p = c(1L, 3L, 1L), rank = 1L, active = 1L),
    list(label = "rectangular", n = 30L, p = c(3L, 4L, 5L), rank = 2L, active = 2L),
    list(label = "wide_training_views", n = 8L, p = c(10L, 3L, 6L), rank = 2L, active = 2L),
    list(label = "zero_spectral_basis", n = 12L, p = c(1L, 3L, 1L), rank = 1L, active = 1L))
  for (rank in c(1L, 2L, 5L)) cases[[length(cases) + 1L]] <- list(
    label = paste0("user_rank_", rank), n = 30L, p = rep(15L, 3L), rank = rank, active = 5L)
  native_rows <- list()
  for (i in seq_along(cases)) {
    one <- cases[[i]]
    sim <- egcar::egcar_simulate(n = one$n, p_list = one$p, rank = one$rank,
      active_per_view = one$active,
      toeplitz_rho = seq(0.5, 0.9, length.out = length(one$p)),
      signal = 0.8, seed = 220L + i)
    if (one$label == "zero_spectral_basis") sim$views[[1L]][, ] <- 0
    prepared <- egcar::egcar_prepare(sim$views)
    z <- e$egcar_get_context(prepared$prep)
    for (group in c(FALSE, TRUE)) {
      family <- if (group) "l21" else "l11"
      s <- e$egcar_initial_state(prepared$prep, z, NULL, group)
      before_context <- serialize(z, NULL); before_state <- serialize(s, NULL)
      controls <- e$egcar_controls(0.005, 0.7, 12L, 0, 0, TRUE,
                                  10, 2, 10L, 5L, group)
      got <- native(z, s, controls, group, FALSE)
      assert(identical(got$matrix_api, 2L),
        "Native matrix API mismatch: reinstall patched egcar and restart R/workers.")
      state_check(got$state, z, group, "cpp")
      reference <- e$egcar_solve_R(z, s, controls, group, FALSE)
      state_check(reference$state, z, group, "r")
      fields <- if (group) c("C", "Gk", "Gl", "Vk", "Vl") else c("C", "Z", "H")
      maximum <- 0
      for (nm in fields) for (j in seq_along(z$edge_k)) {
        a <- got$state[[nm]][[j]]; b <- reference$state[[nm]][[j]]
        matrix_check(a, nrow(b), ncol(b), paste(one$label, family, nm, j))
        maximum <- max(maximum, agree(a, b, paste(one$label, family, nm, j)))
      }
      for (nm in c("primal", "dual", "eps_primal", "eps_dual", "mu"))
        agree(got[[nm]], reference[[nm]], paste(one$label, family, nm))
      assert(identical(got$iterations, reference$iterations) &&
             identical(got$converged, reference$converged),
             paste(one$label, family, "iteration/convergence mismatch."))
      matrix_check(got$history, nrow(reference$history), 7L, "native history")
      agree(got$history, reference$history, paste(one$label, family, "history"))
      assert(identical(serialize(z, NULL), before_context) &&
             identical(serialize(s, NULL), before_state),
             paste(one$label, family, "mutated an input object."))
      native_rows[[length(native_rows) + 1L]] <- data.frame(
        case = one$label, n = one$n, rank = one$rank, penalty = family,
        max_absolute_state_error = maximum, passed = TRUE)
    }
  }
  say("Direct native matrix checks passed: ", length(native_rows), " case/family combinations.")

  fit_rows <- cv_rows <- list()
  run_public_checks <- function() {
    for (rank in c(1L, 2L, 5L)) {
      sim <- egcar::egcar_simulate(n = 30L, p_list = rep(15L, 3L), rank = rank,
        active_per_view = 5L, signal = 0.8, seed = 500L + rank)
      prepared <- egcar::egcar_prepare(sim$views)
      shared <- if (check_cv) egcar::egcar_cv_data(sim$views, nfolds = 5L,
                                                  seed = 600L + rank) else NULL
      for (family in c("l11", "l21")) {
        label <- paste0("r=", rank, ", ", family)
        say("Checking ", label, ": fits, warm start, path, rate",
            if (check_cv) paste0(" and ", workers, "-worker CV") else "")
        fitted <- egcar::egcar_fit(prepared, rank, family, 0.005, ctl)
        before_fit <- serialize(fitted$solver, NULL)
        warm <- egcar::egcar_fit(prepared, rank, family, 0.001, ctl, init = fitted)
        assert(identical(serialize(fitted$solver, NULL), before_fit),
               paste0(label, ": warm start mutated the preceding fit."))
        path <- egcar::egcar_path(prepared, rank, family, c(0.001, 0.005), ctl)
        rate <- egcar::egcar_rate(prepared, rank, family, multiplier = 1, control = ctl)
        fits <- list(direct = fitted, warm = warm, path_1 = path$fits[[1L]],
                     path_2 = path$fits[[2L]], rate = rate)
        for (stage in names(fits)) {
          f <- fits[[stage]]
          # A rate penalty may legitimately select too few rows. That is not a
          # matrix-interface defect; retain and report invalid_loading explicitly.
          check_fit(f, prepared, rank, family, paste(label, stage),
                    require_loading = stage != "rate")
          fit_rows[[length(fit_rows) + 1L]] <<- data.frame(
            rank = rank, signal = 0.8, penalty = family, stage = stage,
            status = f$status, valid_loading = !is.null(f$L), passed = TRUE)
        }
        if (check_cv) {
          serial <- egcar::egcar_cv(shared, rank, family, c(0.001, 0.005),
                                    workers = 1L, control = ctl)
          check_cv_result(serial, paste(label, "serial"))
          check_fit(serial$fit, shared$full, rank, family, paste(label, "serial CV"))
          score_error <- 0
          if (workers > 1L) {
            parallel <- egcar::egcar_cv(shared, rank, family, c(0.001, 0.005),
                                       workers = workers, control = ctl)
            check_cv_result(parallel, paste(label, "parallel"))
            check_fit(parallel$fit, shared$full, rank, family, paste(label, "parallel CV"))
            assert(identical(serial$fold_id, parallel$fold_id) &&
                   identical(serial$cv_table$lambda, parallel$cv_table$lambda) &&
                   identical(serial$lambda, parallel$lambda),
                   paste0(label, ": serial/parallel folds, grids or selected penalty differ."))
            score_error <- agree(serial$cv_fold_table$score,
              parallel$cv_fold_table$score, paste(label, "CV scores"), tolerance = 1e-7)
          }
          cv_rows[[length(cv_rows) + 1L]] <<- data.frame(
            rank = rank, signal = 0.8, penalty = family, folds = 5L,
            workers_checked = workers, candidates = 2L,
            max_absolute_serial_parallel_score_error = score_error, passed = TRUE)
        }
      }
    }
    invisible(NULL)
  }
  if (check_cv && workers > 1L) scope_workers(workers, run_public_checks()) else run_public_checks()
  out <- list(native = do.call(rbind, native_rows), fits = do.call(rbind, fit_rows),
    cv = if (length(cv_rows)) do.call(rbind, cv_rows) else data.frame(),
    package_version = as.character(utils::packageVersion("egcar")),
    package_path = find.package("egcar"), session_info = utils::sessionInfo())
  say("Matrix-interface regression checks passed. Nonconvergence and legitimate ",
      "invalid rate loadings are reported separately; this is not a performance experiment.")
  invisible(out)
}

# Comparison infrastructure and optional dependency diagnostics.

#' Report installed optional comparison dependencies
#'
#' @description Check installed package versions and the expected nonexported PMA helper names without installing packages.
#' @return A data frame of packages, installation status, versions, expected-helper presence and purpose. Helper presence is not a guarantee of numerical compatibility.
#' @rdname benchmark_dependencies
#' @examples
#' benchmark_dependencies()
#' @export
benchmark_dependencies <- function() {
  packages <- c("Rcpp", "RcppArmadillo", "future", "future.apply", "RSpectra",
                "RhpcBLASctl", "RGCCA", "PMA", "ggplot2")
  purpose <- c("Native interface", "Native matrix algebra", "Parallel CV", "Parallel CV",
    "Optional partial loading eigensolve", "Optional scoped BLAS limits",
    "RGCCA and SGCCA", "MultiCCA", "Experiment plots")
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  version <- vapply(seq_along(packages), function(i) {
    if (installed[[i]]) as.character(utils::packageVersion(packages[[i]])) else NA_character_
  }, character(1))
  compatible <- installed
  k <- match("PMA", packages)
  if (installed[[k]]) compatible[[k]] <- all(vapply(c("MultiCCA", "BinarySearch", "soft", "l2n"),
    exists, logical(1), envir = asNamespace("PMA"), inherits = FALSE))
  data.frame(package = packages, installed = installed, version = version,
             expected_helpers_present = compatible, purpose = purpose)
}

# Optional comparison solvers retain the supplied common-loss wrappers.
.egcar_external_cv <- function(x, rank, method, fold_id, nfolds, seed, workers,
                               control, benchmarks, parameters, call) {
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  .egcar_scalar(workers, "workers", 1, integer = TRUE)
  control <- .egcar_as_control(control)
  benchmarks <- .egcar_as_benchmark_control(benchmarks)
  pkg <- switch(method, SGCA = NULL, RGCCA = "RGCCA", SGCCA = "RGCCA", MultiCCA = "PMA")
  if (!is.null(pkg) && !requireNamespace(pkg, quietly = TRUE))
    stop(method, " requires optional package ", pkg, ".", call. = FALSE)
  data <- .egcar_cv_input(x, fold_id, nfolds, seed)
  if (rank > data$p) stop("rank exceeds the total number of variables.")
  workers <- as.integer(min(workers, data$nfolds))
  e <- .egcar_engine(control, benchmarks, workers)
  # Tied bound grids depend on the actual feature dimensions, not a simulation setting.
  if (method == "SGCA" && is.null(parameters$k_grid))
    parameters$k_grid <- sort(unique(c(5L, 10L, 15L, 20L, 30L, data$p)))
  if (method == "SGCCA" && is.null(parameters$sparsity_grid))
    parameters$sparsity_grid <- seq(max(1 / sqrt(data$p_list)), 1, length.out = 10L)
  if (method == "MultiCCA" && is.null(parameters$penalty_grid))
    parameters$penalty_grid <- seq(1, min(sqrt(data$p_list)), length.out = 10L)
  f <- switch(method, SGCA = e$sgca_common_cv, RGCCA = e$rgcca_common_cv,
              SGCCA = e$sgcca_common_cv, MultiCCA = e$multicca_common_cv)
  .egcar_with_seed(seed, .egcar_with_threads(control$blas_threads,
    .egcar_with_workers(workers, {
      out <- do.call(f, c(list(full_views = data$views, full_prep = data$full$prep,
        fold_objects = data$folds, rank = as.integer(rank), seed = seed,
        parallel_folds = workers > 1L), parameters))
      out$method <- method; out$rank <- as.integer(rank)
      out$fold_id <- data$fold_id; out$fold_labels <- data$fold_labels; out$workers <- workers
      out$means <- data$full$means; out$view_names <- data$full$view_names
      out$feature_names <- data$full$feature_names
      out$p_list <- data$p_list; out$p <- data$p; out$n <- data$n
      out$preparation_time <- data$preparation_time
      out$total_time <- out$fit_time + out$tuning_time
      out$loading_time <- 0 # Already included in the external refit timer.
      out$mean_loss <- if (is.null(out$best)) Inf else out$best$mean_loss[[1L]]
      out$call <- call; out$control <- control; out$benchmarks <- benchmarks
      class(out) <- "egcar_cv"
      out
    })))
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

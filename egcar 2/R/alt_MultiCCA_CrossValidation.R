# MultiCCA comparison: common-loss CV and cached-Gram PMA implementation.

#' @rdname comparison_cv
#' @export
multicca_cv <- function(x, rank = 1L, penalty_grid = NULL, fold_id = NULL,
                        nfolds = 5L, seed = 1L, workers = 1L,
                        control = egcar_control(), benchmarks = benchmark_control()) {
  .egcar_external_cv(x, rank, "MultiCCA", fold_id, nfolds, seed, workers, control,
                     benchmarks, list(penalty_grid = penalty_grid), match.call())
}

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

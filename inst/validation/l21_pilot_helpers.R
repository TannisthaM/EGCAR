# Independent calibration helpers. No supplied experiment metrics are read.
# This constructs the same latent-factor distribution as egcar_simulate(),
# omitting population-wide eigenvectors that are unnecessary for fitting/scoring.
l21_pilot_population <- function(p_per_block, rank, signal, seed) {
  ns <- asNamespace("egcar")
  set.seed(seed)
  pp <- rep.int(as.integer(p_per_block), 3L)
  S <- lapply(c(0.5, 0.7, 0.9), function(rho) stats::toeplitz(rho^(0:(p_per_block - 1L))))
  U <- lapply(seq_len(3L), function(k) {
    G <- matrix(0, pp[[k]], rank)
    G[seq_len(5L), ] <- matrix(stats::rnorm(5L * rank), 5L, rank)
    get("metric_normalize", ns)(G, S[[k]])
  })
  B <- lapply(seq_len(3L), function(k) sqrt(signal) * S[[k]] %*% U[[k]])
  Psi <- lapply(seq_len(3L), function(k) get("symmetrize", ns)(
    S[[k]] - signal * S[[k]] %*% U[[k]] %*% t(U[[k]]) %*% S[[k]]))
  list(K = 3L, p_list = pp, rank = rank, signal = signal,
       Sigma_kk = S, U = U, B = B, Psi = Psi)
}

l21_pilot_sample <- function(pop, n, seed) {
  x <- get("simulate_views", asNamespace("egcar"))(pop, n, seed)
  names(x) <- paste0("view", seq_along(x))
  for (k in seq_along(x)) colnames(x[[k]]) <- paste0("V", seq_len(ncol(x[[k]])))
  x
}

l21_pilot_zero_bound <- function(prepared) {
  prep <- prepared$prep
  accum <- lapply(prep$p_list, numeric)
  for (j in seq_len(nrow(prep$edge_table))) {
    e <- prep$edge_table[j, ]
    S <- prep$S_kl[[e$key]]
    accum[[e$k]] <- accum[[e$k]] + rowSums(S * S)
    accum[[e$l]] <- accum[[e$l]] + colSums(S * S)
  }
  0.5 * sqrt(max(unlist(accum, use.names = FALSE)))
}

l21_pilot_metric_error <- function(L, pop) {
  if (is.null(L)) return(NA_real_)
  ii <- get("make_block_indices", asNamespace("egcar"))(pop$p_list)
  H <- matrix(0, pop$rank, pop$rank)
  D <- H
  for (k in seq_len(pop$K)) {
    A <- L[ii[[k]], , drop = FALSE]
    H <- H + crossprod(A, pop$Sigma_kk[[k]] %*% A)
    D <- D + crossprod(A, pop$Sigma_kk[[k]] %*% pop$U[[k]]) / sqrt(pop$K)
  }
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE)
  if (min(ev$values) <= 1e-12 * max(ev$values)) return(NA_real_)
  normalized <- sweep(crossprod(ev$vectors, D), 1L, sqrt(ev$values), "/")
  sqrt(max(0, pop$rank - sum(pmin(1, svd(normalized, nu = 0, nv = 0)$d)^2)))
}

l21_pilot_design <- function(stage = c("pilot", "confirmation")) {
  stage <- match.arg(stage)
  ranks <- if (stage == "pilot") c(1L, 3L, 5L) else seq_len(5L)
  d <- expand.grid(p_per_block = c(10L, 100L, 500L, 1000L),
                   rank = ranks, signal = c(0.3, 0.5, 0.8), KEEP.OUT.ATTRS = FALSE)
  d <- d[order(d$p_per_block, d$signal, d$rank), ]
  rownames(d) <- NULL
  d$scenario <- seq_len(nrow(d))
  d$n_train <- 150L
  d$n_validation <- 600L
  d$population_seed <- (if (stage == "pilot") 100214000L else 200214000L) + 1000L * d$scenario
  d$training_seed <- d$population_seed + 77L
  d$validation_seed <- d$population_seed + 177L
  d$stage <- stage
  d
}

l21_pilot_evaluate <- function(row, multipliers, control) {
  pop <- l21_pilot_population(row$p_per_block, row$rank, row$signal, row$population_seed)
  train <- l21_pilot_sample(pop, row$n_train, row$training_seed)
  validation <- l21_pilot_sample(pop, row$n_validation, row$validation_seed)
  prep <- egcar::egcar_prepare(train)
  base <- sqrt((max(prep$p - prep$p_list) + log(prep$p)) / prep$n)
  bound <- l21_pilot_zero_bound(prep)
  out <- vector("list", length(multipliers))
  for (j in seq_along(multipliers)) {
    cg <- multipliers[[j]]
    lambda <- cg * base
    start <- proc.time()[[3L]]
    certified <- lambda > bound
    if (certified) {
      one <- list(status = "invalid_loading", converged = TRUE, iterations = 0L,
                  L = NULL, selected = integer(), error = "Certified unique zero optimum")
      score <- NA_real_
    } else {
      one <- tryCatch(egcar::egcar_rate(prep, row$rank, penalty = "l21",
                                       multiplier = cg, control = control),
                      error = function(e) list(status = "error", converged = FALSE,
                        iterations = NA_integer_, L = NULL, selected = integer(),
                        error = conditionMessage(e)))
      score <- if (identical(one$status, "ok")) egcar::egcar_score(one, validation) else NA_real_
    }
    # Truth is only used after the multiplier has been frozen, for confirmation.
    error <- if (row$stage == "confirmation") l21_pilot_metric_error(one$L, pop) else NA_real_
    out[[j]] <- cbind(row, data.frame(multiplier = cg, lambda = lambda,
      zero_bound = bound, zero_certified = certified, status = one$status,
      converged = one$converged, iterations = one$iterations,
      selected_rows = length(one$selected), validation_score = score,
      normalized_score = score / row$rank, subspace_sigma0 = error,
      elapsed = proc.time()[[3L]] - start,
      error_message = if (is.null(one$error)) NA_character_ else one$error,
      stringsAsFactors = FALSE))
    cat(sprintf("%s %02d: p/block=%d r=%d signal=%.1f c=%.2f status=%s rows=%d score/r=%s\n",
      row$stage, row$scenario, row$p_per_block, row$rank, row$signal, cg,
      one$status, length(one$selected), format(score / row$rank, digits = 5)))
    flush.console()
    rm(one); invisible(gc(FALSE))
  }
  do.call(rbind, out)
}

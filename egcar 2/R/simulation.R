# The unchanged positive-semidefinite multiview simulation model.

#' Simulation model, subspace distance and oracle diagnostics
#'
#' @description Reproduce the supplied positive-semidefinite multiview latent-factor model and its two truth-informed oracle benchmarks.
#' @param n Number of observations.
#' @param p_list Positive view dimensions.
#' @param rank Latent rank or requested loading rank.
#' @param active_per_view Active variable count per view, restricted by each view dimension.
#' @param toeplitz_rho One Toeplitz correlation or one per view, all with absolute value below one.
#' @param signal Shared signal strictly between zero and one.
#' @param seed Nonnegative integer seed; caller RNG state is restored.
#' @param A,B Finite matrices with equal row counts; their column spaces are compared.
#' @param x Raw views, prepared data or shared CV data for the support oracle.
#' @param active_local List of known true local column supports. Must not be used as information in ordinary CV.
#' @param ridge Small support-oracle linear-solve ridge.
#' @param control An \code{egcar_control} object or named subset.
#' @param population Population component returned by \code{egcar_simulate}.
#' @details The model, covariance divisor, normalization and oracle definitions are
#' inherited from the supplied script. Oracle1 runs the original zero-penalty consensus
#' ADMM on exact population covariance blocks with both penalty coefficients equal
#' to zero; this diagnostic is intentionally not substituted with the accelerated
#' group-only splitting. Oracle2 solves on empirical covariance blocks restricted
#' to known true supports. Neither oracle performs cross-validation.
#' \code{egcar_distance} is the Frobenius sine-theta distance after Euclidean
#' orthonormalization, not an aligned coefficient difference. It returns NA for
#' numerically deficient rank. Premultiply both bases by the population metric
#' square root for the supplied Sigma0-subspace diagnostic.
#' @return Simulation returns \code{views}, \code{population}, \code{truth} and support information. Distance returns a scalar. Oracles return fitted objects.
#' @rdname egcar_simulation
#' @examples
#' sim <- egcar_simulate(n = 30)
#' oracle <- egcar_oracle_support(sim$views, sim$active_local)
#' egcar_distance(oracle$L, sim$truth)
#' @export
egcar_simulate <- function(n = 120L, p_list = c(15L, 15L, 15L), rank = 1L,
                           active_per_view = 5L, toeplitz_rho = c(0.5, 0.7, 0.9),
                           signal = 0.8, seed = 1L) {
  .egcar_scalar(n, "n", 3, integer = TRUE)
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  .egcar_scalar(active_per_view, "active_per_view", 1, integer = TRUE)
  if (!is.numeric(p_list) || length(p_list) < 2L || any(!is.finite(p_list)) ||
      any(p_list < 1 | p_list != floor(p_list)) || any(p_list > .Machine$integer.max))
    stop("p_list must contain at least two positive integers.")
  if (rank > min(active_per_view, p_list)) stop("rank must not exceed the active size of any view.")
  if (length(toeplitz_rho) == 1L) toeplitz_rho <- rep(toeplitz_rho, length(p_list))
  if (!is.numeric(toeplitz_rho) || length(toeplitz_rho) != length(p_list) ||
      any(!is.finite(toeplitz_rho)) || any(abs(toeplitz_rho) >= 1))
    stop("toeplitz_rho must contain one value or one per view, all with absolute value < 1.")
  .egcar_scalar(signal, "signal", strict = TRUE)
  if (signal >= 1) stop("signal must be less than one.")
  .egcar_scalar(seed, "seed", 0, integer = TRUE)
  .egcar_with_seed(seed, {
    population <- make_population(as.integer(p_list), as.integer(rank), as.integer(active_per_view),
                                  toeplitz_rho, signal, as.integer(seed))
    next_seed <- as.integer((as.double(seed) + 77) %% 2147483646)
    views <- simulate_views(population, as.integer(n), next_seed)
    names(views) <- paste0("view", seq_along(views))
    views <- .egcar_views(views)
    list(views = views, population = population, truth = population$Lstar,
         active_local = population$active_local, seed = seed)
  })
}

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

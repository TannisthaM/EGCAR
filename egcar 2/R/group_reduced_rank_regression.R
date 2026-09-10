# View-wise l21-only EGCAR and incident-edge row-group operations.

fit_l21_admm <- function(prep, lambda_g, mu_g = 1, max_iter = 2000L,
    abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
    scale_factor = 2, adapt_every = 10L, group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10, init = NULL, keep_history = FALSE,
    verbose = FALSE, check_every = 1L) {
  if (EGCAR_BACKEND == "reference") return(fit_l21_admm_reference(
    prep, lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu, balance_ratio,
    scale_factor, adapt_every, group_zero_tol, entry_zero_tol, init,
    keep_history, verbose, check_every))
  ctl <- egcar_controls(lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu,
    balance_ratio, scale_factor, adapt_every, check_every, keep_history)
  a <- egcar_run_solver(prep, ctl, init, TRUE, verbose)
  s <- a$state; z <- a$context; names(s$C) <- z$keys
  G <- egcar_view_copies(z, s$Gk, s$Gl)
  V <- egcar_view_copies(z, s$Vk, s$Vl)
  active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
  C_hat <- lapply(s$C, function(A) { A[abs(A) < entry_zero_tol] <- 0; A })
  for (e in seq_along(z$edge_k)) {
    C_hat[[e]][!active_rows[[z$edge_k[[e]]]], ] <- 0
    C_hat[[e]][, !active_rows[[z$edge_l[[e]]]]] <- 0
  }
  history <- if (keep_history) {
    h <- as.data.frame(a$history)
    names(h) <- c("iter", "objective", "primal", "dual", "eps_primal", "eps_dual", "mu_g")
    h$iter <- as.integer(h$iter); h
  } else NULL
  list(C_hat = C_hat, C = s$C, G = G, V = V, active_rows = active_rows,
    converged = a$converged, iterations = a$iterations, primal_residual = a$primal,
    dual_residual = a$dual, eps_primal = a$eps_primal, eps_dual = a$eps_dual,
    rho_e = 0, lambda_g = lambda_g, mu_g = a$mu, history = history)
}

egcar_view_copies <- function(z, left, right) {
  p <- sum(z$p_list)
  out <- lapply(z$p_list, function(pk) matrix(0, pk, p - pk))
  for (e in seq_along(z$edge_k)) {
    out[[z$edge_k[[e]]]][, z$cols_k[[e]]] <- left[[e]]
    out[[z$edge_l[[e]]]][, z$cols_l[[e]]] <- t(right[[e]])
  }
  out
}

egcar_group_norms <- function(z, left, right = left) {
  ans <- lapply(z$p_list, numeric)
  for (e in seq_along(z$edge_k)) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    ans[[k]] <- ans[[k]] + rowSums(left[[e]] * left[[e]])
    ans[[l]] <- ans[[l]] + colSums(right[[e]] * right[[e]])
  }
  lapply(ans, sqrt)
}

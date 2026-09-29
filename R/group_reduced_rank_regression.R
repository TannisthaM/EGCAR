# View-wise l21-only EGCAR and incident-edge row-group operations.

fit_l21_admm <- function(prep, lambda_g, mu_g = 1, max_iter = 2000L,
    abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
    scale_factor = 2, adapt_every = 10L, group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10, init = NULL, keep_history = FALSE,
    verbose = FALSE, check_every = 1L, keep_state = TRUE, compact_state = FALSE) {
  if (EGCAR_BACKEND == "reference") {
    if (all(c("Hk", "Hl", "a") %in% names(init)))
      init <- egcar_expand_group_state(egcar_get_context(prep), init)
    if (all(c("Gk", "Gl", "Vk", "Vl") %in% names(init))) {
      z <- egcar_get_context(prep)
      init$G <- egcar_view_copies(z, init$Gk, init$Gl)
      init$V <- egcar_view_copies(z, init$Vk, init$Vl)
    }
    out <- fit_l21_admm_reference(
      prep, lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu, balance_ratio,
      scale_factor, adapt_every, group_zero_tol, entry_zero_tol, init,
      keep_history, verbose, check_every)
    if (!keep_state) out[c("C", "G", "V")] <- NULL
    return(out)
  }
  ctl <- egcar_controls(lambda_g, mu_g, max_iter, abs_tol, rel_tol, adaptive_mu,
    balance_ratio, scale_factor, adapt_every, check_every, keep_history)
  a <- egcar_run_solver(prep, ctl, init, TRUE, verbose)
  s <- a$state; z <- a$context; names(s$C) <- z$keys
  # Compute active rows directly from endpoint group copies. The larger
  # view-wide G/V matrices are only materialized when warm-start state is needed.
  active_rows <- lapply(Map(`*`, egcar_group_norms(z, s$Hk, s$Hl), s$a),
                        function(x) x > group_zero_tol)
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
  out <- list(C_hat = C_hat, active_rows = active_rows,
    converged = a$converged, iterations = a$iterations, primal_residual = a$primal,
    dual_residual = a$dual, eps_primal = a$eps_primal, eps_dual = a$eps_dual,
    rho_e = 0, lambda_g = lambda_g, mu_g = a$mu, history = history)
  if (keep_state) {
    out$C <- s$C
    if (compact_state) {
      fields <- c("Hk", "Hl", "a")
      out[fields] <- s[fields]
    } else {
      legacy <- egcar_expand_group_state(z, s)
      out$G <- egcar_view_copies(z, legacy$Gk, legacy$Gl)
      out$V <- egcar_view_copies(z, legacy$Vk, legacy$Vl)
    }
  }
  out
}

# Explicit compatibility conversion; ordinary fits and CV retain Hk/Hl/a.
egcar_expand_group_state <- function(z, s) {
  if (is.null(s[["a"]])) return(s)
  E <- seq_along(z$edge_k)
  list(C = s$C,
    Gk = lapply(E, function(e) s$Hk[[e]] * s$a[[z$edge_k[[e]]]]),
    Gl = lapply(E, function(e) sweep(s$Hl[[e]], 2L, s$a[[z$edge_l[[e]]]], "*")),
    Vk = lapply(E, function(e) s$Hk[[e]] * (1 - s$a[[z$edge_k[[e]]]])),
    Vl = lapply(E, function(e) sweep(s$Hl[[e]], 2L, 1 - s$a[[z$edge_l[[e]]]], "*")))
}

# Coerce one edge block to its mathematically known p_k x p_l shape.
# A dimensionless numeric vector is accepted only when its length is exactly
# p_k*p_l. R and Armadillo are both column-major, so matrix() restores the
# legacy nested-Rcpp output without reordering values. Any other shape fails.
egcar_edge_matrix <- function(z, A, e, label = "edge block") {
  k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
  nr <- as.integer(z$p_list[[k]]); nc <- as.integer(z$p_list[[l]])
  expected <- c(nr, nc)
  if (is.matrix(A) && is.numeric(A) && identical(dim(A), expected)) return(A)
  if (is.numeric(A) && is.null(dim(A)) && length(A) == nr * nc) {
    return(matrix(as.numeric(A), nrow = nr, ncol = nc))
  }
  observed <- if (is.null(dim(A))) {
    paste0(typeof(A), " object of length ", length(A), " with no dim attribute")
  } else {
    paste0(typeof(A), " object with dimensions ", paste(dim(A), collapse = " x "))
  }
  stop("EGCAR matrix interface: ", label, " for edge ", k, "_", l,
       " must be a numeric ", nr, " x ", nc, " matrix; got ", observed, ".",
       call. = FALSE)
}

egcar_view_copies <- function(z, left, right) {
  p <- sum(z$p_list)
  out <- lapply(z$p_list, function(pk) matrix(0, pk, p - pk))
  for (e in seq_along(z$edge_k)) {
    L <- egcar_edge_matrix(z, left[[e]], e, "left endpoint copy")
    R <- egcar_edge_matrix(z, right[[e]], e, "right endpoint copy")
    out[[z$edge_k[[e]]]][, z$cols_k[[e]]] <- L
    out[[z$edge_l[[e]]]][, z$cols_l[[e]]] <- t(R)
  }
  out
}

egcar_group_norms <- function(z, left, right = left) {
  if (EGCAR_BACKEND == "cpp") {
    # Preserve the compatibility shape checks without copying valid matrices.
    left <- lapply(seq_along(z$edge_k), function(e) egcar_edge_matrix(z, left[[e]], e))
    right <- lapply(seq_along(z$edge_k), function(e) egcar_edge_matrix(z, right[[e]], e))
    return(egcar_native_group_norms(left, right, z$edge_k, z$edge_l, z$p_list))
  }
  ans <- lapply(z$p_list, numeric)
  for (e in seq_along(z$edge_k)) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    L <- egcar_edge_matrix(z, left[[e]], e, "left group block")
    R <- egcar_edge_matrix(z, right[[e]], e, "right group block")
    ans[[k]] <- ans[[k]] + rowSums(L * L)
    ans[[l]] <- ans[[l]] + colSums(R * R)
  }
  lapply(ans, sqrt)
}

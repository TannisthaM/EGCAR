# Unchanged dense-reference solvers and original zero-penalty consensus oracle solver.

cache_problem_matrices_reference <- function(prep) {
  prep$indices <- make_block_indices(prep$p_list)
  prep$Sigma0 <- block_diag(prep$S_kk)
  prep$eig_products <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    outer(prep$eig[[k]]$values, prep$eig[[l]]$values, "*")
  })
  prep$edge_cols_k <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
    k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
    prep$layout[[k]]$cols[[as.character(l)]]
  })
  prep$edge_cols_l <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
    k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
    prep$layout[[l]]$cols[[as.character(k)]]
  })
  prep$loading_factor_cache <- new.env(parent = emptyenv())
  prep
}

prepare_problem_reference <- function(centered_views) {
  K <- length(centered_views)
  n <- nrow(centered_views[[1L]])
  p_list <- vapply(centered_views, ncol, integer(1L))
  if (any(vapply(centered_views, nrow, integer(1L)) != n)) {
    stop("All views must have the same number of rows.")
  }

  edge_table <- make_edge_table(p_list)
  layout <- make_incidence_layout(p_list)
  S_kk <- lapply(centered_views, function(X) crossprod(X) / n)
  S_kl <- setNames(vector("list", nrow(edge_table)), edge_table$key)

  for (e in seq_len(nrow(edge_table))) {
    k <- edge_table$k[[e]]
    l <- edge_table$l[[e]]
    S_kl[[edge_table$key[[e]]]] <-
      crossprod(centered_views[[k]], centered_views[[l]]) / n
  }

  eig <- lapply(S_kk, function(S) {
    ee <- eigen(symmetrize(S), symmetric = TRUE)
    ee$values <- pmax(ee$values, 0)
    ee
  })

  cache_problem_matrices(list(
    n = n,
    K = K,
    p_list = p_list,
    p = sum(p_list),
    edge_table = edge_table,
    layout = layout,
    S_kk = S_kk,
    S_kl = S_kl,
    eig = eig,
    q = sum(edge_table$p_k * edge_table$p_l)
  ))
}

fit_l11_admm_reference <- function(
    prep,
    rho_e,
    mu = 1,
    max_iter = 2000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    adaptive_mu = TRUE,
    balance_ratio = 10,
    scale_factor = 2,
    adapt_every = 10L,
    entry_zero_tol = 1e-10,
    init = NULL,
    verbose = FALSE,
    check_every = 1L) {

  C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
  Z <- if (!is.null(init$Z)) init$Z else empty_edge_list(prep$edge_table)
  H <- if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)
  converged <- FALSE
  r_primal <- Inf
  r_dual <- Inf
  eps_primal <- NA_real_
  eps_dual <- NA_real_

  if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
  edge_k <- prep$edge_table$k
  edge_l <- prep$edge_table$l
  edge_keys <- prep$edge_table$key
  n_edges <- length(edge_keys)
  Q <- lapply(prep$eig, `[[`, "vectors")
  denominator_cache <- NULL
  previous_shift <- NA_real_

  for (iter in seq_len(max_iter)) {
    current_shift <- mu
    if (!identical(current_shift, previous_shift)) {
      denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
      previous_shift <- current_shift
    }
    Z_old <- Z
    C_new <- C  # all edges are overwritten; no zero matrices allocated

    for (e in seq_len(n_edges)) {
      k <- edge_k[[e]]; l <- edge_l[[e]]
      B <- prep$S_kl[[e]] + mu * (Z[[e]] - H[[e]])
      C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
      C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
    }

    Z_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
      key <- prep$edge_table$key[[e]]
      soft_threshold(C_new[[key]] + H[[key]], rho_e / mu)
    }), prep$edge_table$key)

    H_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
      key <- prep$edge_table$key[[e]]
      H[[key]] + C_new[[key]] - Z_new[[key]]
    }), prep$edge_table$key)

    C <- C_new
    Z <- Z_new
    H <- H_new

    # See fit_oracle_consensus_admm for rationale: skip residual/convergence/
    # adaptive-mu bookkeeping on non-check iterations to cut overhead
    # across the very large number of solves used in cross-validation.
    do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

    if (do_check) {
      r_primal <- sqrt(sum(vapply(prep$edge_table$key, function(key) {
        sum((C_new[[key]] - Z_new[[key]])^2)
      }, numeric(1L))))
      r_dual <- mu * sqrt(sum(vapply(prep$edge_table$key, function(key) {
        sum((Z_new[[key]] - Z_old[[key]])^2)
      }, numeric(1L))))

      norm_C <- sqrt(sum(vapply(C_new, function(A) sum(A^2), numeric(1L))))
      norm_Z <- sqrt(sum(vapply(Z_new, function(A) sum(A^2), numeric(1L))))
      eps_primal <- sqrt(prep$q) * abs_tol + rel_tol * max(norm_C, norm_Z)
      norm_dual <- mu * sqrt(sum(vapply(H_new, function(A) sum(A^2), numeric(1L))))
      eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * norm_dual

      if (r_primal <= eps_primal && r_dual <= eps_dual) {
        converged <- TRUE
        break
      }

      if (adaptive_mu && iter %% adapt_every == 0L) {
        if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
          mu <- mu * scale_factor
          H <- lapply(H, function(A) A / scale_factor)
        } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
          mu <- mu / scale_factor
          H <- lapply(H, function(A) A * scale_factor)
        }
      }

      if (verbose && (iter == 1L || iter %% 100L == 0L)) {
        cat(sprintf("  pairwise iter=%d primal=%.3e dual=%.3e\n", iter, r_primal, r_dual))
      }
    }
  }

  C_hat <- lapply(Z, function(A) {
    A[abs(A) < entry_zero_tol] <- 0
    A
  })

  list(
    C_hat = C_hat,
    C = C,
    Z = Z,
    H = H,
    converged = converged,
    iterations = iter,
    primal_residual = r_primal,
    dual_residual = r_dual,
    eps_primal = eps_primal,
    eps_dual = eps_dual,
    rho_e = rho_e,
    lambda_g = 0,
    mu_z = mu
  )
}

fit_l21_admm_reference <- function(
    prep,
    lambda_g,
    mu_g = 1,
    max_iter = 2000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    adaptive_mu = TRUE,
    balance_ratio = 10,
    scale_factor = 2,
    adapt_every = 10L,
    group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10,
    init = NULL,
    keep_history = FALSE,
    verbose = FALSE,
    check_every = 1L) {

  if (length(lambda_g) != 1L || !is.finite(lambda_g) || lambda_g < 0) {
    stop("lambda_g must be a finite, nonnegative scalar.")
  }
  if (length(mu_g) != 1L || !is.finite(mu_g) || mu_g <= 0) {
    stop("mu_g must be a finite, positive scalar.")
  }

  C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)

  M0 <- assemble_all_M(C, prep$layout)
  G <- if (!is.null(init$G)) init$G else lapply(M0, function(A) A)
  V <- if (!is.null(init$V)) init$V else lapply(M0, function(A) matrix(0, nrow(A), ncol(A)))

  history <- if (keep_history) {
    data.frame(
      iter = integer(0), objective = numeric(0),
      primal = numeric(0), dual = numeric(0),
      eps_primal = numeric(0), eps_dual = numeric(0),
      mu_g = numeric(0)
    )
  } else NULL

  converged <- FALSE
  r_primal <- Inf
  r_dual <- Inf
  eps_primal <- NA_real_
  eps_dual <- NA_real_

  if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
  edge_k <- prep$edge_table$k
  edge_l <- prep$edge_table$l
  edge_keys <- prep$edge_table$key
  n_edges <- length(edge_keys)
  Q <- lapply(prep$eig, `[[`, "vectors")
  denominator_cache <- NULL
  previous_shift <- NA_real_

  for (iter in seq_len(max_iter)) {
    # Each edge occurs once at each endpoint: M* M = 2 I.
    current_shift <- 2 * mu_g
    if (!identical(current_shift, previous_shift)) {
      denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
      previous_shift <- current_shift
    }
    G_old <- G
    C_new <- C  # all edges are overwritten; no zero matrices allocated
    # Each endpoint difference is formed once, not once for every incident edge.
    GV <- Map(function(A, B) A - B, G, V)

    # C updates retain exactly the same Sylvester solution in cached eigenbases.
    for (e in seq_len(n_edges)) {
      k <- edge_k[[e]]; l <- edge_l[[e]]
      target_k <- GV[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
      target_l <- t(GV[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
      B <- prep$S_kl[[e]] + mu_g * (target_k + target_l)
      C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
      C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
    }

    M_new <- assemble_all_M(C_new, prep$layout)

    G_new <- lapply(seq_len(prep$K), function(k) {
      row_group_threshold(M_new[[k]] + V[[k]], lambda_g / mu_g)
    })

    V_new <- lapply(seq_len(prep$K), function(k) {
      V[[k]] + M_new[[k]] - G_new[[k]]
    })

    C <- C_new
    G <- G_new
    V <- V_new

    # Residuals, tolerances, adaptive-mu, and convergence are only
    # recomputed every `check_every` iterations (always including the
    # first and last iteration). The ADMM updates above still run every
    # iteration; only this bookkeeping is skipped on non-check iterations,
    # which removes most of the R-level overhead across the very large
    # number of solves performed during cross-validation.
    do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

    if (do_check) {
      r_primal_sq <- 0
      for (k in seq_len(prep$K)) {
        r_primal_sq <- r_primal_sq + sum((M_new[[k]] - G_new[[k]])^2)
      }
      r_primal <- sqrt(r_primal_sq)

      dG <- Map(function(A, B) A - B, G_new, G_old)
      r_dual_sq <- 0
      for (e in seq_len(nrow(prep$edge_table))) {
        k <- prep$edge_table$k[[e]]
        l <- prep$edge_table$l[[e]]
        key <- prep$edge_table$key[[e]]
        dGk <- dG[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        dGl <- t(dG[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        S_edge <- mu_g * (dGk + dGl)
        r_dual_sq <- r_dual_sq + sum(S_edge^2)
      }
      r_dual <- sqrt(r_dual_sq)

      norm_Ax_sq <- sum(vapply(M_new, function(A) sum(A^2), numeric(1L)))
      norm_Bz_sq <- sum(vapply(G_new, function(A) sum(A^2), numeric(1L)))
      d_primal <- sum(prep$p_list * (prep$p - prep$p_list))  # 2 * prep$q
      eps_primal <- sqrt(d_primal) * abs_tol +
        rel_tol * max(sqrt(norm_Ax_sq), sqrt(norm_Bz_sq))

      dual_adj_sq <- 0
      for (e in seq_len(nrow(prep$edge_table))) {
        k <- prep$edge_table$k[[e]]
        l <- prep$edge_table$l[[e]]
        key <- prep$edge_table$key[[e]]
        Vk <- V_new[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        Vl <- t(V_new[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        dual_adj <- mu_g * (Vk + Vl)
        dual_adj_sq <- dual_adj_sq + sum(dual_adj^2)
      }
      eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * sqrt(dual_adj_sq)

      if (keep_history) {
        objective_value <- operator_objective(prep, C_new, rho_e = 0, lambda_g = lambda_g)
        history <- rbind(
          history,
          data.frame(
            iter = iter,
            objective = objective_value,
            primal = r_primal,
            dual = r_dual,
            eps_primal = eps_primal,
            eps_dual = eps_dual,
            mu_g = mu_g
          )
        )
      }

      if (!is.finite(r_primal) || !is.finite(r_dual)) {
        warning("ADMM produced a non-finite residual.")
        break
      }

      if (r_primal <= eps_primal && r_dual <= eps_dual) {
        converged <- TRUE
        break
      }

      if (adaptive_mu && iter %% adapt_every == 0L) {
        if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
          mu_g <- mu_g * scale_factor
          V <- lapply(V, function(A) A / scale_factor)
        } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
          mu_g <- mu_g / scale_factor
          V <- lapply(V, function(A) A * scale_factor)
        }
      }

      if (verbose && (iter == 1L || iter %% 100L == 0L)) {
        cat(sprintf(
          "  iter=%d  primal=%.3e (%.3e)  dual=%.3e (%.3e)\n",
          iter, r_primal, eps_primal, r_dual, eps_dual
        ))
      }
    }
  }

  # Return C, masking rows/columns selected out by the endpoint groups.
  # The existing entry_zero_tol is retained ONLY as a numerical cutoff;
  # there is no entrywise penalty, Z variable, or soft-thresholding step.
  active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
  C_hat <- lapply(C, function(A) {
    A[abs(A) < entry_zero_tol] <- 0
    A
  })
  for (e in seq_len(nrow(prep$edge_table))) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    key <- prep$edge_table$key[[e]]
    C_hat[[key]][!active_rows[[k]], ] <- 0
    C_hat[[key]][, !active_rows[[l]]] <- 0
  }

  list(
    C_hat = C_hat,
    C = C,
    G = G,
    V = V,
    active_rows = active_rows,
    converged = converged,
    iterations = iter,
    primal_residual = r_primal,
    dual_residual = r_dual,
    eps_primal = eps_primal,
    eps_dual = eps_dual,
    rho_e = 0,  # compatibility/reporting only; no entrywise penalty
    lambda_g = lambda_g,
    mu_g = mu_g,
    history = history
  )
}

fit_oracle_consensus_admm <- function(
    prep,
    mu_z = 1,
    mu_g = 1,
    max_iter = 2000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    adaptive_mu = TRUE,
    balance_ratio = 10,
    scale_factor = 2,
    adapt_every = 10L,
    group_zero_tol = 1e-8,
    entry_zero_tol = 1e-10,
    init = NULL,
    keep_history = FALSE,
    verbose = FALSE,
    check_every = 1L) {

  # Oracle1 only: both coefficients are fixed at zero, never tuned.
  rho_e <- lambda_g <- 0
  if (mu_z <= 0 || mu_g <= 0) stop("ADMM augmentation parameters must be positive.")

  C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
  Z <- if (!is.null(init$Z)) init$Z else lapply(C, function(A) A)
  H <- if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)

  M0 <- assemble_all_M(C, prep$layout)
  G <- if (!is.null(init$G)) init$G else lapply(M0, function(A) A)
  V <- if (!is.null(init$V)) init$V else lapply(M0, function(A) matrix(0, nrow(A), ncol(A)))

  history <- if (keep_history) {
    data.frame(
      iter = integer(0), objective = numeric(0),
      primal = numeric(0), dual = numeric(0),
      eps_primal = numeric(0), eps_dual = numeric(0),
      mu_z = numeric(0), mu_g = numeric(0)
    )
  } else NULL

  converged <- FALSE
  r_primal <- Inf
  r_dual <- Inf
  eps_primal <- NA_real_
  eps_dual <- NA_real_

  if (is.null(prep$eig_products)) prep <- cache_problem_matrices(prep)
  edge_k <- prep$edge_table$k
  edge_l <- prep$edge_table$l
  edge_keys <- prep$edge_table$key
  n_edges <- length(edge_keys)
  Q <- lapply(prep$eig, `[[`, "vectors")
  denominator_cache <- NULL
  previous_shift <- NA_real_

  for (iter in seq_len(max_iter)) {
    current_shift <- mu_z + 2 * mu_g
    if (!identical(current_shift, previous_shift)) {
      denominator_cache <- lapply(prep$eig_products, function(D) D + current_shift)
      previous_shift <- current_shift
    }
    Z_old <- Z
    G_old <- G
    C_new <- C  # all edges are overwritten; no zero matrices allocated
    # Each endpoint difference is formed once, not once for every incident edge.
    GV <- Map(function(A, B) A - B, G, V)

    # C updates retain exactly the same Sylvester solution in cached eigenbases.
    for (e in seq_len(n_edges)) {
      k <- edge_k[[e]]; l <- edge_l[[e]]
      target_k <- GV[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
      target_l <- t(GV[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
      B <- prep$S_kl[[e]] + mu_z * (Z[[e]] - H[[e]]) + mu_g * (target_k + target_l)
      C_tilde <- (crossprod(Q[[k]], B) %*% Q[[l]]) / denominator_cache[[e]]
      C_new[[e]] <- tcrossprod(Q[[k]] %*% C_tilde, Q[[l]])
    }

    M_new <- assemble_all_M(C_new, prep$layout)

    Z_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
      key <- prep$edge_table$key[[e]]
      soft_threshold(C_new[[key]] + H[[key]], rho_e / mu_z)
    }), prep$edge_table$key)

    G_new <- lapply(seq_len(prep$K), function(k) {
      row_group_threshold(M_new[[k]] + V[[k]], lambda_g / mu_g)
    })

    H_new <- setNames(lapply(seq_len(nrow(prep$edge_table)), function(e) {
      key <- prep$edge_table$key[[e]]
      H[[key]] + C_new[[key]] - Z_new[[key]]
    }), prep$edge_table$key)

    V_new <- lapply(seq_len(prep$K), function(k) {
      V[[k]] + M_new[[k]] - G_new[[k]]
    })

    C <- C_new
    Z <- Z_new
    G <- G_new
    H <- H_new
    V <- V_new

    # Residuals, tolerances, adaptive-mu, and convergence are only
    # recomputed every `check_every` iterations (always including the
    # first and last iteration). The ADMM updates above still run every
    # iteration; only this bookkeeping is skipped on non-check iterations,
    # which removes most of the R-level overhead across the very large
    # number of solves performed during cross-validation.
    do_check <- (iter %% check_every == 0L) || (iter == 1L) || (iter == max_iter)

    if (do_check) {
      r_primal_sq <- 0
      for (key in prep$edge_table$key) {
        r_primal_sq <- r_primal_sq + sum((C_new[[key]] - Z_new[[key]])^2)
      }
      for (k in seq_len(prep$K)) {
        r_primal_sq <- r_primal_sq + sum((M_new[[k]] - G_new[[k]])^2)
      }
      r_primal <- sqrt(r_primal_sq)

      dG <- Map(function(A, B) A - B, G_new, G_old)
      r_dual_sq <- 0
      for (e in seq_len(nrow(prep$edge_table))) {
        k <- prep$edge_table$k[[e]]
        l <- prep$edge_table$l[[e]]
        key <- prep$edge_table$key[[e]]
        dGk <- dG[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        dGl <- t(dG[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        S_edge <- mu_z * (Z_new[[key]] - Z_old[[key]]) + mu_g * (dGk + dGl)
        r_dual_sq <- r_dual_sq + sum(S_edge^2)
      }
      r_dual <- sqrt(r_dual_sq)

      norm_Ax_sq <- sum(vapply(C_new, function(A) sum(A^2), numeric(1L))) +
        sum(vapply(M_new, function(A) sum(A^2), numeric(1L)))
      norm_Bz_sq <- sum(vapply(Z_new, function(A) sum(A^2), numeric(1L))) +
        sum(vapply(G_new, function(A) sum(A^2), numeric(1L)))
      d_primal <- prep$q + sum(prep$p_list * (prep$p - prep$p_list))
      eps_primal <- sqrt(d_primal) * abs_tol +
        rel_tol * max(sqrt(norm_Ax_sq), sqrt(norm_Bz_sq))

      dual_adj_sq <- 0
      for (e in seq_len(nrow(prep$edge_table))) {
        k <- prep$edge_table$k[[e]]
        l <- prep$edge_table$l[[e]]
        key <- prep$edge_table$key[[e]]
        Vk <- V_new[[k]][, prep$edge_cols_k[[e]], drop = FALSE]
        Vl <- t(V_new[[l]][, prep$edge_cols_l[[e]], drop = FALSE])
        dual_adj <- mu_z * H_new[[key]] + mu_g * (Vk + Vl)
        dual_adj_sq <- dual_adj_sq + sum(dual_adj^2)
      }
      eps_dual <- sqrt(prep$q) * abs_tol + rel_tol * sqrt(dual_adj_sq)

      if (keep_history) {
        objective_value <- operator_objective(prep, C_new, rho_e, lambda_g)
        history <- rbind(
          history,
          data.frame(
            iter = iter,
            objective = objective_value,
            primal = r_primal,
            dual = r_dual,
            eps_primal = eps_primal,
            eps_dual = eps_dual,
            mu_z = mu_z,
            mu_g = mu_g
          )
        )
      }

      if (!is.finite(r_primal) || !is.finite(r_dual)) {
        warning("ADMM produced a non-finite residual.")
        break
      }

      if (r_primal <= eps_primal && r_dual <= eps_dual) {
        converged <- TRUE
        break
      }

      if (adaptive_mu && iter %% adapt_every == 0L) {
        if (r_primal > balance_ratio * max(r_dual, .Machine$double.eps)) {
          mu_z <- mu_z * scale_factor
          mu_g <- mu_g * scale_factor
          H <- lapply(H, function(A) A / scale_factor)
          V <- lapply(V, function(A) A / scale_factor)
        } else if (r_dual > balance_ratio * max(r_primal, .Machine$double.eps)) {
          mu_z <- mu_z / scale_factor
          mu_g <- mu_g / scale_factor
          H <- lapply(H, function(A) A * scale_factor)
          V <- lapply(V, function(A) A * scale_factor)
        }
      }

      if (verbose && (iter == 1L || iter %% 100L == 0L)) {
        cat(sprintf(
          "  iter=%d  primal=%.3e (%.3e)  dual=%.3e (%.3e)\n",
          iter, r_primal, eps_primal, r_dual, eps_dual
        ))
      }
    }
  }

  # Enforce both proximal sparsity patterns in the reported estimate.  This
  # changes only the final tolerance-level disagreement between consensus copies.
  active_rows <- lapply(G, function(A) row_l2(A) > group_zero_tol)
  C_hat <- lapply(Z, function(A) {
    A[abs(A) < entry_zero_tol] <- 0
    A
  })
  for (e in seq_len(nrow(prep$edge_table))) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    key <- prep$edge_table$key[[e]]
    C_hat[[key]][!active_rows[[k]], ] <- 0
    C_hat[[key]][, !active_rows[[l]]] <- 0
  }

  list(
    C_hat = C_hat,
    C = C,
    Z = Z,
    G = G,
    H = H,
    V = V,
    active_rows = active_rows,
    converged = converged,
    iterations = iter,
    primal_residual = r_primal,
    dual_residual = r_dual,
    eps_primal = eps_primal,
    eps_dual = eps_dual,
    rho_e = rho_e,
    lambda_g = lambda_g,
    mu_z = mu_z,
    mu_g = mu_g,
    history = history
  )
}

loading_from_operator <- function(
    prep,
    C,
    rank,
    row_threshold = 1e-4,
    covariance_ridge = 1e-4,
    require_positive = TRUE,
    positive_tol = 1e-10) {

  C_full <- assemble_full_C(C, prep$p_list)
  selected <- which(row_l2(C_full) > row_threshold)
  if (length(selected) < rank) {
    return(list(valid = FALSE, reason = "fewer selected rows than rank"))
  }

  factors <- loading_metric_factors(prep, selected, covariance_ridge)
  S_half <- factors$half
  S_inv_half <- factors$inv_half
  R_sel <- symmetrize(
    S_half %*% C_full[selected, selected, drop = FALSE] %*% S_half
  )
  ee <- eigen(R_sel, symmetric = TRUE)
  if (length(ee$values) < rank) {
    return(list(valid = FALSE, reason = "operator dimension below rank"))
  }
  if (require_positive && ee$values[[rank]] <= positive_tol) {
    return(list(valid = FALSE, reason = "fewer than rank positive eigenvalues"))
  }

  U <- ee$vectors[, seq_len(rank), drop = FALSE]
  L_sel <- S_inv_half %*% U
  L <- matrix(0, prep$p, rank)
  L[selected, ] <- L_sel

  list(
    valid = TRUE,
    L = L,
    U = U,
    selected = selected,
    eigenvalues = ee$values[seq_len(rank)],
    generalized_eigenvalues = 1 + ee$values[seq_len(rank)],
    C_full = C_full
  )
}

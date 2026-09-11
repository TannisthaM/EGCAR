# Entrywise-only EGCAR and shared accelerated ADMM dispatch.

fit_l11_admm <- function(prep, rho_e, mu = 1, max_iter = 2000L,
    abs_tol = 1e-5, rel_tol = 1e-4, adaptive_mu = TRUE, balance_ratio = 10,
    scale_factor = 2, adapt_every = 10L, entry_zero_tol = 1e-10,
    init = NULL, verbose = FALSE, check_every = 1L, keep_state = TRUE) {
  if (EGCAR_BACKEND == "reference") {
    out <- fit_l11_admm_reference(
      prep, rho_e, mu, max_iter, abs_tol, rel_tol, adaptive_mu, balance_ratio,
      scale_factor, adapt_every, entry_zero_tol, init, verbose, check_every)
    if (!keep_state) out[c("C", "Z", "H")] <- NULL
    return(out)
  }
  ctl <- egcar_controls(rho_e, mu, max_iter, abs_tol, rel_tol, adaptive_mu,
    balance_ratio, scale_factor, adapt_every, check_every, FALSE)
  a <- egcar_run_solver(prep, ctl, init, FALSE, verbose)
  s <- a$state; keys <- a$context$keys
  for (nm in c("C", "Z", "H")) names(s[[nm]]) <- keys
  C_hat <- lapply(s$Z, function(A) { A[abs(A) < entry_zero_tol] <- 0; A })
  out <- list(C_hat = C_hat,
    converged = a$converged, iterations = a$iterations, primal_residual = a$primal,
    dual_residual = a$dual, eps_primal = a$eps_primal, eps_dual = a$eps_dual,
    rho_e = rho_e, lambda_g = 0, mu_z = a$mu)
  if (keep_state) out[c("C", "Z", "H")] <- list(s$C, s$Z, s$H)
  out
}

egcar_solve_R <- function(z, s, ctl, group, verbose = FALSE) {
  E <- seq_along(z$edge_k); K <- seq_along(z$p_list)
  mu <- ctl$mu; last_shift <- NA_real_; den <- NULL
  converged <- FALSE; rp <- rd <- Inf; ep <- ed <- NA_real_
  history <- if (ctl$history && group) matrix(NA_real_, ctl$max_iter, 7L) else NULL
  h <- 0L
  Ct <- Wk <- Wl <- vector("list", length(E))
  for (it in seq_len(ctl$max_iter)) {
    check <- it == 1L || it == ctl$max_iter || it %% ctl$check_every == 0L
    shift <- if (group) 2 * mu else mu
    if (!identical(shift, last_shift)) {
      den <- lapply(z$D, function(D) D + shift); last_shift <- shift
    }
    rp2 <- rd2 <- nc2 <- nz2 <- ny2 <- 0
    for (e in E) {
      T <- if (group) (s$Gk[[e]] - s$Vk[[e]] + s$Gl[[e]] - s$Vl[[e]]) / 2 else
        s$Z[[e]] - s$H[[e]]
      Pt <- egcar_project(z, e, T)
      Ct[[e]] <- (z$St[[e]] + shift * Pt) / den[[e]]
      s$C[[e]] <- if (z$full[[e]]) egcar_lift(z, e, Ct[[e]]) else
        T + z$remainder[[e]] / shift + egcar_lift(z, e, Ct[[e]] - Pt)
      if (!group) {
        W <- s$C[[e]] + s$H[[e]]
        Zn <- soft_threshold(W, ctl$penalty / mu)
        Hn <- W - Zn
        if (check) {
          rp2 <- rp2 + sum((s$C[[e]] - Zn)^2)
          rd2 <- rd2 + sum((Zn - s$Z[[e]])^2)
          nc2 <- nc2 + sum(s$C[[e]]^2); nz2 <- nz2 + sum(Zn^2)
          ny2 <- ny2 + sum(Hn^2)
        }
        s$Z[[e]] <- Zn; s$H[[e]] <- Hn
      } else {
        Wk[[e]] <- s$C[[e]] + s$Vk[[e]]
        Wl[[e]] <- s$C[[e]] + s$Vl[[e]]
      }
    }
    if (group) {
      norms <- egcar_group_norms(z, Wk, Wl)
      threshold <- ctl$penalty / mu
      mult <- if (threshold <= 0) lapply(z$p_list, function(pk) rep.int(1, pk)) else
        lapply(norms, function(nr) pmax(0, 1 - threshold / pmax(nr, .Machine$double.eps)))
      for (e in E) {
        k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
        Gkn <- Wk[[e]] * mult[[k]]
        # sweep scales columns; no transpose of the large edge is required.
        Gln <- sweep(Wl[[e]], 2L, mult[[l]], "*")
        Vkn <- Wk[[e]] - Gkn; Vln <- Wl[[e]] - Gln
        if (check) {
          rp2 <- rp2 + sum((s$C[[e]] - Gkn)^2) + sum((s$C[[e]] - Gln)^2)
          rd2 <- rd2 + sum((Gkn - s$Gk[[e]] + Gln - s$Gl[[e]])^2)
          nc2 <- nc2 + 2 * sum(s$C[[e]]^2)
          nz2 <- nz2 + sum(Gkn^2) + sum(Gln^2)
          ny2 <- ny2 + sum((Vkn + Vln)^2)
        }
        s$Gk[[e]] <- Gkn; s$Gl[[e]] <- Gln
        s$Vk[[e]] <- Vkn; s$Vl[[e]] <- Vln
      }
    }
    if (check) {
      rp <- sqrt(rp2); rd <- mu * sqrt(rd2)
      ep <- sqrt(if (group) 2 * z$q else z$q) * ctl$abs_tol +
        ctl$rel_tol * max(sqrt(nc2), sqrt(nz2))
      ed <- sqrt(z$q) * ctl$abs_tol + ctl$rel_tol * mu * sqrt(ny2)
      if (ctl$history && group) {
        objective <- 0
        for (e in E) objective <- objective +
          0.5 * sum(z$D[[e]] * Ct[[e]]^2) - sum(z$S[[e]] * s$C[[e]])
        objective <- objective + ctl$penalty * sum(unlist(egcar_group_norms(z, s$C), use.names = FALSE))
        h <- h + 1L; history[h, ] <- c(it, objective, rp, rd, ep, ed, mu)
      }
      if (!is.finite(rp) || !is.finite(rd)) { warning("ADMM produced a non-finite residual."); break }
      if (rp <= ep && rd <= ed) { converged <- TRUE; break }
      if (ctl$adaptive && it %% ctl$adapt_every == 0L) {
        factor <- if (rp > ctl$balance_ratio * max(rd, .Machine$double.eps)) ctl$scale_factor else
          if (rd > ctl$balance_ratio * max(rp, .Machine$double.eps)) 1 / ctl$scale_factor else 1
        if (factor != 1) {
          mu <- mu * factor
          if (group) {
            s$Vk <- lapply(s$Vk, function(A) A / factor)
            s$Vl <- lapply(s$Vl, function(A) A / factor)
          } else s$H <- lapply(s$H, function(A) A / factor)
        }
      }
      if (verbose && (it == 1L || it %% 100L == 0L)) cat(sprintf(
        "  iter=%d primal=%.3e (%.3e) dual=%.3e (%.3e)\n", it, rp, ep, rd, ed))
    }
  }
  list(state = s, converged = converged, iterations = it, primal = rp, dual = rd,
       eps_primal = ep, eps_dual = ed, mu = mu,
       history = if (is.null(history)) matrix(numeric(), 0L, 7L) else history[seq_len(h), , drop = FALSE])
}

egcar_controls <- function(penalty, mu, max_iter, abs_tol, rel_tol, adaptive_mu,
                           balance_ratio, scale_factor, adapt_every, check_every, keep_history) {
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(penalty) || penalty < 0 || !scalar(mu) || mu <= 0)
    stop("EGCAR requires a finite nonnegative penalty and a finite positive ADMM parameter.")
  if (!scalar(abs_tol) || abs_tol < 0 || !scalar(rel_tol) || rel_tol < 0)
    stop("EGCAR tolerances must be finite and nonnegative.")
  for (x in list(max_iter, adapt_every, check_every))
    if (!scalar(x) || x < 1 || x != floor(x)) stop("EGCAR iteration counts must be positive integers.")
  if (!scalar(balance_ratio) || balance_ratio <= 0 || !scalar(scale_factor) || scale_factor <= 1)
    stop("Invalid residual-balancing parameters.")
  list(penalty = penalty, mu = mu, max_iter = as.integer(max_iter), abs_tol = abs_tol,
    rel_tol = rel_tol, adaptive = isTRUE(adaptive_mu), balance_ratio = balance_ratio,
    scale_factor = scale_factor, adapt_every = as.integer(adapt_every),
    check_every = as.integer(check_every), history = isTRUE(keep_history))
}

egcar_run_solver <- function(prep, ctl, init, group, verbose = FALSE) {
  z <- egcar_get_context(prep)
  s <- egcar_initial_state(prep, z, init, group)
  raw <- if (EGCAR_BACKEND == "cpp") {
    egcar_native_solve(z, s, ctl, group, verbose)
  } else {
    egcar_solve_R(z, s, ctl, group, verbose)
  }
  # Defensive compatibility layer: restore only a missing dim attribute using
  # the exact edge dimensions known from the prepared problem. The current
  # native API already returns matrices explicitly, so this is normally a no-op.
  raw$state <- .egcar_normalize_solver_state(raw$state, z, group)
  if (EGCAR_BACKEND == "cpp" && !is.null(raw$matrix_api) && !identical(raw$matrix_api, 2L))
    stop("EGCAR native matrix API mismatch. Reinstall the patched egcar source ",
         "package and restart R, including CV workers.", call. = FALSE)
  .egcar_check_solver_state(raw$state, z, group, EGCAR_BACKEND)
  raw$context <- z
  raw
}

.egcar_solve_one <- function(prep, penalty, lambda, control, e, init = NULL,
                             cv = FALSE, history = control$keep_history, keep_state = TRUE) {
  init <- .egcar_validate_init(init, prep, penalty)
  args <- list(prep = prep, max_iter = if (cv) control$max_iter_cv else control$max_iter,
    abs_tol = control$abs_tol, rel_tol = control$rel_tol,
    adaptive_mu = control$adaptive_mu, balance_ratio = control$balance_ratio,
    scale_factor = control$scale_factor, adapt_every = control$adapt_every,
    entry_zero_tol = control$entry_zero_tol, init = init,
    verbose = control$verbose,
    check_every = if (cv) control$check_every_cv else control$check_every)
  if (penalty == "l11") {
    args$rho_e <- lambda; args$mu <- control$mu; args$keep_state <- keep_state
    do.call(e$fit_l11_admm, args)
  } else {
    args$lambda_g <- lambda; args$mu_g <- control$mu
    args$group_zero_tol <- control$group_zero_tol; args$keep_history <- history
    args$keep_state <- keep_state
    do.call(e$fit_l21_admm, args)
  }
}

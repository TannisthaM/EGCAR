# Truth-informed population and support oracles; never used for ordinary tuning.

#' @rdname egcar_simulation
#' @export
egcar_oracle_support <- function(x, active_local, rank = 1L, ridge = 1e-8,
                                 control = egcar_control()) {
  prepared <- egcar_prepare(x)
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  .egcar_scalar(ridge, "ridge")
  if (rank > prepared$p) stop("rank exceeds p.")
  if (!is.list(active_local) || length(active_local) != length(prepared$p_list))
    stop("active_local must contain one support vector per view.")
  for (k in seq_along(active_local)) {
    a <- active_local[[k]]
    if (!is.numeric(a) || !length(a) || any(!is.finite(a)) ||
        any(a != floor(a) | a < 1 | a > prepared$p_list[[k]]) || anyDuplicated(a))
      stop("Each oracle support must be a nonempty vector of distinct valid column indices.")
  }
  control <- .egcar_as_control(control)
  e <- .egcar_engine(control)
  .egcar_with_threads(control$blas_threads, {
    t0 <- proc.time()[[3L]]
    C <- e$fit_oracle_support(prepared$prep, active_local, ridge)
    fit_time <- proc.time()[[3L]] - t0
    raw <- list(C_hat = C, C = C, converged = TRUE, iterations = 1L)
    t0 <- proc.time()[[3L]]
    loading <- .egcar_loading(prepared$prep, raw, as.integer(rank), control, e)
    out <- .egcar_fit_object(raw, loading, prepared, as.integer(rank), "oracle-support", 0,
                             control, fit_time, proc.time()[[3L]] - t0, match.call())
    out$active_local <- active_local
    out
  })
}

#' @rdname egcar_simulation
#' @export
egcar_oracle_population <- function(population, rank = population$rank,
                                    control = egcar_control()) {
  required <- c("p_list", "Sigma_kk", "Sigma_kl", "rank")
  if (!is.list(population) || !all(required %in% names(population)))
    stop("population must be the population component returned by egcar_simulate().")
  .egcar_scalar(rank, "rank", 1, integer = TRUE)
  control <- .egcar_as_control(control)
  e <- .egcar_engine(control)
  # Preserve Oracle1 exactly: original consensus splitting with BOTH coefficients 0.
  # This is not the accelerated l21-only estimator and is never cross-validated.
  .egcar_with_seed(NULL, .egcar_with_threads(control$blas_threads, {
    prep <- e$prepare_population_problem(population)
    if (rank > prep$p) stop("rank exceeds p.")
    raw <- e$fit_oracle_population(population, as.integer(rank), prep,
                                    keep_history = control$keep_history)
    prepared <- list(p_list = prep$p_list, p = prep$p, n = Inf,
      means = lapply(prep$p_list, function(p) rep(0, p)),
      view_names = paste0("view", seq_along(prep$p_list)),
      feature_names = lapply(prep$p_list, function(p) paste0("V", seq_len(p))),
      preparation_time = 0)
    if (is.null(raw$raw_fit)) stop(raw$error)
    out <- .egcar_fit_object(raw$raw_fit, raw$loading, prepared, as.integer(rank),
                             "oracle-population", 0, control, raw$time, 0, match.call())
    out$oracle_definition <- "Population covariances, original zero-penalty consensus splitting, both coefficients zero."
    out
  }))
}

prepare_population_problem <- function(population) {
  p_list <- population$p_list
  edge_table <- make_edge_table(p_list)
  eig <- lapply(population$Sigma_kk, function(S) {
    ee <- eigen(symmetrize(S), symmetric = TRUE)
    ee$values <- pmax(ee$values, 0)
    ee
  })
  cache_problem_matrices(list(
    n = Inf,
    K = length(p_list),
    p_list = p_list,
    p = sum(p_list),
    edge_table = edge_table,
    layout = make_incidence_layout(p_list),
    S_kk = population$Sigma_kk,
    S_kl = population$Sigma_kl,
    eig = eig,
    q = sum(edge_table$p_k * edge_table$p_l)
  ))
}

fit_zero_penalty_oracle_admm <- function(
    prep,
    rank,
    covariance_ridge,
    max_iter,
    source_label,
    keep_history = TRUE,
    retain_raw_fit = TRUE) {

  start <- proc.time()[[3L]]

  out <- tryCatch({
    fit <- fit_oracle_consensus_admm(
      prep = prep,
      mu_z = MU_Z,
      mu_g = MU_G,
      max_iter = max_iter,
      abs_tol = ABS_TOL,
      rel_tol = REL_TOL,
      adaptive_mu = ADAPTIVE_MU,
      group_zero_tol = GROUP_ZERO_TOL,
      entry_zero_tol = ENTRY_ZERO_TOL,
      keep_history = keep_history
    )

    loading <- loading_from_operator(
      prep = prep,
      C = fit$C_hat,
      rank = rank,
      row_threshold = ROW_THRESHOLD,
      covariance_ridge = covariance_ridge,
      require_positive = TRUE
    )
    if (!retain_raw_fit && is.list(loading)) loading$C_full <- NULL

    messages <- character(0)
    status <- "ok"
    if (!isTRUE(fit$converged)) {
      status <- "not_converged"
      messages <- c(messages, sprintf(
        "%s zero-penalty oracle ADMM did not satisfy the stopping rule in %d iterations.",
        source_label, fit$iterations
      ))
    }
    if (!isTRUE(loading$valid)) {
      if (identical(status, "ok")) status <- "invalid_loading"
      messages <- c(messages, paste0(
        "Loading extraction failed: ", loading$reason %||% "unknown reason", "."
      ))
    }

    list(
      C = fit$C_hat,
      loading = loading,
      prep = if (retain_raw_fit) prep else NULL,
      raw_fit = if (retain_raw_fit) fit else NULL,
      status = status,
      error = if (length(messages) == 0L) NA_character_ else paste(messages, collapse = " "),
      converged = fit$converged,
      iterations = fit$iterations,
      rho_e = 0,
      lambda_g = 0
    )
  }, error = function(e) {
    list(
      C = NULL,
      loading = NULL,
      prep = prep,
      raw_fit = NULL,
      status = "failed",
      error = conditionMessage(e),
      converged = FALSE,
      iterations = NA_integer_,
      rho_e = 0,
      lambda_g = 0
    )
  })

  out$time <- proc.time()[[3L]] - start
  out
}

fit_oracle_population <- function(
    population,
    rank,
    prep = NULL,
    keep_history = TRUE,
    retain_raw_fit = TRUE) {

  if (is.null(prep)) prep <- prepare_population_problem(population)
  fit_zero_penalty_oracle_admm(
    prep = prep,
    rank = rank,
    covariance_ridge = 0,
    max_iter = ORACLE1_MAX_ITER,
    source_label = "Population",
    keep_history = keep_history,
    retain_raw_fit = retain_raw_fit
  )
}

fit_oracle_support <- function(prep, active_local, ridge = 1e-8) {
  C <- empty_edge_list(prep$edge_table)
  for (e in seq_len(nrow(prep$edge_table))) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    key <- prep$edge_table$key[[e]]
    sk <- active_local[[k]]
    sl <- active_local[[l]]
    A <- prep$S_kk[[k]][sk, sk, drop = FALSE] + ridge * diag(length(sk))
    B <- prep$S_kk[[l]][sl, sl, drop = FALSE] + ridge * diag(length(sl))
    R <- prep$S_kl[[key]][sk, sl, drop = FALSE]
    C_sub <- t(solve(B, t(solve(A, R))))
    C[[key]][sk, sl] <- C_sub
  }
  C
}

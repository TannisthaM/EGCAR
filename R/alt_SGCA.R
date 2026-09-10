# SGCA comparison: common-loss CV, cached initializer and penalized TGD.

#' Common-loss CV for SGCA, RGCCA, SGCCA and MultiCCA
#'
#' @description Optional comparison wrappers retaining the supplied script's fitting algorithms, common validation score and fold-level parallelism.
#' @param x Raw views or the same \code{egcar_cv_data} object used for EGCAR.
#' @param rank Target loading rank.
#' @param k_grid SGCA total row-support budget across all views. NULL uses valid members of 5, 10, 15, 20, 30 and p.
#' @param rho_grid SGCA direct initializer coefficients; not rate multipliers.
#' @param lambda_grid SGCA strictly positive penalized-gradient coefficients.
#' @param tau_grid RGCCA tied tau grid in [0,1].
#' @param sparsity_grid SGCCA tied sparsity bounds between the maximum of 1/sqrt(p_k) and one. NULL uses a dimension-valid grid.
#' @param penalty_grid MultiCCA tied L1 bounds between one and the minimum of sqrt(p_k). NULL uses a dimension-valid grid.
#' @param fold_id,nfolds,seed,workers Common-fold and parallel controls as in \code{egcar_cv}.
#' @param control EGCAR infrastructure controls, including scoped BLAS threads. EGCAR ADMM controls do not change comparison solvers.
#' @param benchmarks A \code{benchmark_control} object or named argument list.
#' @details SGCA uses the bundled initializer
#' \code{sgca_init_fixed}, from the TannisthaM/SGCA repository. Its initializer is
#' cached by fold and direct rho. The subsequent penalized thresholded-gradient
#' implementation is the one in the supplied script, including backtracking and
#' no per-iteration metric normalization; the supplied local-experiment CV and TGD routines are retained.
#' It tunes k, rho and lambda jointly.
#'
#' RGCCA and SGCCA use the public \code{RGCCA::rgcca} solver with variable/block
#' scaling disabled. They use \code{astar} loadings mapping original centered data
#' to deflated components. MultiCCA uses the supplied cached-Gram implementation
#' and the installed PMA threshold helpers, or \code{PMA::MultiCCA} when selected.
#' External block signs are resolved using training observations only.
#'
#' All wrappers minimize exactly the same held-out loss as EGCAR. The wrappers do
#' not use the comparison packages' own cross-validation criteria. A candidate
#' must have finite loss on every fold. A returned finite approximate fit does
#' not imply convergence; diagnostic flags are retained. These methods return
#' loading estimates, not an EGCAR cross-view operator.
#'
#' Optional packages are never automatically installed. Run the included
#' \code{00_install_dependencies.R --benchmarks} example explicitly when needed.
#' For external methods loading extraction is included in \code{fit_time};
#' \code{loading_time} is therefore zero, not an additional omitted cost.
#' @return An \code{egcar_cv} object with \code{L}, \code{best}, the external \code{fit_full}, candidate and fold diagnostics, status, means, controls, labels and timings.
#' @rdname comparison_cv
#' @examples
#' \dontrun{
#' sim <- egcar_simulate(n = 60)
#' shared <- egcar_cv_data(sim$views, nfolds = 3)
#' rgcca_cv(shared, tau_grid = c(0.1, 1))
#' sgcca_cv(shared, sparsity_grid = c(0.5, 1))
#' multicca_cv(shared, penalty_grid = c(2, 3))
#' sgca_cv(shared, k_grid = c(10, 15), rho_grid = 0.01, lambda_grid = 0.1)
#' }
#' @export
sgca_cv <- function(x, rank = 1L, k_grid = NULL,
                    rho_grid = c(0, 1e-3, 1e-2, 0.1, 0.5, 1),
                    lambda_grid = 10^seq(-5, 4), fold_id = NULL, nfolds = 5L,
                    seed = 1L, workers = 1L, control = egcar_control(),
                    benchmarks = benchmark_control()) {
  .egcar_external_cv(x, rank, "SGCA", fold_id, nfolds, seed, workers, control,
    benchmarks, list(k_grid = k_grid, rho_grid = rho_grid, lambda_grid = lambda_grid), match.call())
}

sgca_hard_rows <- function(U, k) {
  U <- as.matrix(U)
  if (k < ncol(U) || k > nrow(U)) stop("SGCA requires rank <= k <= p.")
  if (k < nrow(U)) {
    keep <- order(rowSums(U * U), decreasing = TRUE)[seq_len(k)]
    U[setdiff(seq_len(nrow(U)), keep), ] <- 0
  }
  U
}

sgca_metric_normalize <- function(U, B, tol = 1e-10) {
  G <- symmetrize(crossprod(U, B %*% U))
  ev <- eigen(G, symmetric = TRUE)
  d <- ev$values
  if (any(!is.finite(d)) || d[[1L]] <= 0 ||
      min(d) <= tol * d[[1L]]) stop("SGCA has a rank-deficient training Gram matrix.")
  tcrossprod(sweep(U %*% ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
}

sgca_prepare_tgd_start <- function(A, B, init, k) {
  U0 <- sgca_metric_normalize(sgca_hard_rows(init, k), B)
  ev <- eigen(symmetrize(crossprod(U0, A %*% U0)), symmetric = TRUE)
  list(U0 = U0, vectors = ev$vectors, values = ev$values,
       U0_vectors = U0 %*% ev$vectors)
}

sgca_prepare_initializer <- function(A, B, reference, nu = 1) {
  A <- symmetrize(A)
  B <- symmetrize(B)
  h_update <- tryCatch(get("updateH", envir = environment(reference), inherits = TRUE),
                       error = function(e) NULL)
  if (!is.function(h_update)) return(NULL)  # keep compatibility with other versions
  ev <- eigen(B, symmetric = TRUE)
  d <- pmax(ev$values, 0)
  sqB <- tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors)
  tau <- 4 * nu * max(d)^2
  if (!is.finite(tau) || tau <= 0) tau <- 1
  list(A = A, B = B, sqB = sqB, tau = tau, nu = nu,
       A_scaled = (1 / tau) * A, B_scaled = (nu / tau) * B,
       sqB_scaled = (nu / tau) * sqB, update_H = h_update)
}

sgca_init_cached <- function(prepared, rho, K,
                             epsilon = SGCA_INIT_TOL,
                             maxiter = SGCA_MAX_ITER_INIT, trace = FALSE) {
  if (is.null(prepared)) stop("Missing SGCA initializer preparation.")
  z <- prepared
  p <- nrow(z$B)
  H <- Pi <- oldPi <- diag(1, p)
  Gamma <- matrix(0, p, p)
  criteria <- Inf
  iter <- 0L
  threshold <- rho / z$tau
  while (criteria > epsilon && iter < maxiter) {
    fixed_H <- z$sqB_scaled %*% (H - Gamma) %*% z$sqB
    for (j in seq_len(20L)) {
      Pi <- soft_threshold(Pi + z$A_scaled -
        z$B_scaled %*% Pi %*% z$B + fixed_H, threshold)
    }
    H <- z$update_H(z$sqB, Gamma, z$nu, Pi, K)
    Gamma <- Gamma + (z$sqB %*% Pi %*% z$sqB - H) * z$nu
    criteria <- sqrt(sum((Pi - oldPi)^2))
    oldPi <- Pi
    iter <- iter + 1L
    if (trace) cat("iter:", iter, "crit:", criteria, "\n")
  }
  list(Pi = Pi, H = H, Gamma = Gamma, iteration = iter, convergence = criteria)
}

sgca_tgd_penalized <- function(
    A, B, init, rank, k, lambda,
    eta = SGCA_ETA, max_iter = SGCA_MAX_ITER_TGD,
    tol = SGCA_TGD_TOL, max_backtrack = 50L,
    prepared_start = NULL, matrices_prepared = FALSE) {

  if (!is.finite(lambda) || lambda <= 0) stop("SGCA lambda must be strictly positive.")
  if (!is.finite(eta) || eta <= 0 || max_iter < 1L) stop("Invalid SGCA TGD controls.")
  if (!isTRUE(matrices_prepared)) {
    A <- symmetrize(A)
    B <- symmetrize(B)
  }
  if (is.null(prepared_start)) prepared_start <- sgca_prepare_tgd_start(A, B, init, k)
  st <- prepared_start
  W <- tcrossprod(sweep(st$U0_vectors, 2L,
                         sqrt(pmax(st$values + lambda, 0)), "*"), st$vectors)
  lambda_I <- lambda * diag(rank)
  objective <- function(Z, AZ, BZ) {
    -0.5 * sum(Z * AZ) + 0.25 * sum((crossprod(Z, BZ) - lambda_I)^2)
  }
  AW <- A %*% W
  BW <- B %*% W
  value <- objective(W, AW, BW)
  if (!is.finite(value)) stop("Non-finite SGCA initialization objective.")
  converged <- FALSE
  step <- 2 * eta
  total_backtracks <- 0L
  relative_mapping <- Inf

  for (iter in seq_len(max_iter)) {
    grad <- -AW + BW %*% (crossprod(W, BW) - lambda_I)
    if (any(!is.finite(grad))) stop("Non-finite SGCA gradient.")
    trial_step <- min(2 * eta, 1.25 * step)
    accepted <- FALSE
    for (bt in 0:max_backtrack) {
      Wnew <- sgca_hard_rows(W - trial_step * grad, k)
      if (all(is.finite(Wnew))) {
        AWnew <- A %*% Wnew
        BWnew <- B %*% Wnew
        new_value <- objective(Wnew, AWnew, BWnew)
        diff_sq <- sum((Wnew - W)^2)
        # Sufficient descent for the projected/hard-thresholded step.
        slack <- 1e-12 * max(1, abs(value))
        if (is.finite(new_value) &&
            new_value <= value - diff_sq / (4 * trial_step) + slack) {
          accepted <- TRUE
          break
        }
      }
      trial_step <- trial_step / 2
    }
    if (!accepted) stop("SGCA backtracking could not find a finite descent step.")
    total_backtracks <- total_backtracks + bt
    relative_mapping <- sqrt(diff_sq) / (trial_step * max(1, frob(W)))
    relative_change <- sqrt(diff_sq) / max(1, frob(W))
    W <- Wnew
    AW <- AWnew
    BW <- BWnew
    value <- new_value
    step <- trial_step
    # A small step alone is NOT interpreted as convergence.
    if (relative_mapping <= tol && relative_change <= tol) {
      converged <- TRUE
      break
    }
  }
  L <- sgca_metric_normalize(W, B)
  L <- validate_loading_matrix(L, nrow(A), rank)
  list(L = L, converged = converged, iterations = iter,
       objective_scaled = value, final_step = step,
       relative_gradient_mapping = relative_mapping,
       backtracking_steps = total_backtracks,
       lambda = lambda, k = k)
}

get_sgca_initializer <- function() {
  # Bundled initializer and its dependency closure: no optional namespace lookup.
  sgca_init_fixed
}

sgca_common_cv <- function(
    full_views, full_prep, fold_objects, rank,
    k_grid = SGCA_K_GRID, rho_grid = SGCA_RHO_GRID,
    lambda_grid = SGCA_LAMBDA_GRID, seed = 1L,
    parallel_folds = PARALLEL_CV) {

  p <- full_prep$p
  if (any(!is.finite(k_grid)) || any(k_grid != floor(k_grid))) {
    stop("SGCA k_grid must contain finite integers.")
  }
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[k_grid >= rank & k_grid <= p]
  if (!length(k_grid)) stop("No SGCA k satisfies rank <= k <= p.")
  if (!length(rho_grid) || any(!is.finite(rho_grid)) || any(rho_grid < 0)) {
    stop("SGCA rho_grid must contain nonnegative, finite DIRECT coefficients.")
  }
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    stop("SGCA lambda_grid must contain positive, finite coefficients.")
  }
  grid <- expand.grid(
    sgca_k = k_grid, sgca_rho = sort(unique(rho_grid), decreasing = TRUE),
    sgca_lambda = sort(unique(lambda_grid)), KEEP.OUT.ATTRS = FALSE
  )
  # Exact CV ties: smaller k, then larger rho, then smaller lambda.
  grid <- grid[order(grid$sgca_k, -grid$sgca_rho, grid$sgca_lambda), , drop = FALSE]
  rownames(grid) <- NULL

  prepare_context <- function(train_views, prep, final) {
    A <- full_covariance_from_prep(prep)
    B <- prep$Sigma0 %||% block_diag(prep$S_kk)
    diag(B) <- diag(B) + SGCA_RIDGE_B
    initializer <- get_sgca_initializer()
    list(
      prep = prep, A = A, B = B, initializer = initializer,
      initializer_prepared = if (FAST_SGCA_INITIALIZER)
        sgca_prepare_initializer(A, B, initializer) else NULL,
      init_cache = new.env(parent = emptyenv()),
      tgd_start_cache = new.env(parent = emptyenv())
    )
  }

  fit_candidate <- function(context, parameters, final) {
    rho <- parameters$sgca_rho[[1L]]
    k <- parameters$sgca_k[[1L]]
    lambda <- parameters$sgca_lambda[[1L]]
    key <- sprintf("rho_%.17g", rho)
    if (!exists(key, envir = context$init_cache, inherits = FALSE)) {
      initialized <- tryCatch({
        z <- if (!is.null(context$initializer_prepared)) {
          sgca_init_cached(context$initializer_prepared, rho = rho, K = rank,
            epsilon = SGCA_INIT_TOL, maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
        } else {
          context$initializer(A = context$A, B = context$B, rho = rho, K = rank,
            nu = 1, epsilon = SGCA_INIT_TOL,
            maxiter = SGCA_MAX_ITER_INIT, trace = FALSE)
        }
        Pi <- as.matrix(z$Pi)
        if (any(!is.finite(Pi))) stop("Non-finite SGCA initializer.")
        ss <- svd(Pi, nu = rank, nv = 0L)
        if (length(ss$d) < rank || ss$d[[1L]] <= 0 ||
            ss$d[[rank]] <= 1e-10 * ss$d[[1L]]) {
          stop("The SGCA initializer has numerical rank below the requested rank.")
        }
        U <- sweep(ss$u[, seq_len(rank), drop = FALSE], 2L,
                   sqrt(ss$d[seq_len(rank)]), "*")
        list(U = U, raw = z, error = NULL)
      }, error = function(e) list(U = NULL, error = conditionMessage(e)))
      assign(key, initialized, envir = context$init_cache)
    }
    ini <- get(key, envir = context$init_cache, inherits = FALSE)
    if (is.null(ini$U)) stop(ini$error)
    start_key <- paste0(key, "_k", k)
    if (!exists(start_key, envir = context$tgd_start_cache, inherits = FALSE)) {
      st <- tryCatch(list(value = sgca_prepare_tgd_start(context$A, context$B, ini$U, k)),
                     error = function(e) list(error = conditionMessage(e)))
      assign(start_key, st, envir = context$tgd_start_cache)
    }
    st <- get(start_key, envir = context$tgd_start_cache, inherits = FALSE)
    if (is.null(st$value)) stop(st$error)
    tgd <- sgca_tgd_penalized(
      A = context$A, B = context$B, init = ini$U,
      rank = rank, k = k, lambda = lambda,
      prepared_start = st$value, matrices_prepared = TRUE
    )
    init_conv <- if (is.numeric(ini$raw$convergence) && length(ini$raw$convergence) == 1L) {
      is.finite(ini$raw$convergence) && ini$raw$convergence <= SGCA_INIT_TOL
    } else NA
    conv <- if (identical(init_conv, FALSE)) FALSE else tgd$converged
    init_iters <- as.integer(ini$raw$iteration %||% NA_integer_)
    list(
      L = tgd$L,
      fit = list(initializer = ini$raw, tgd = tgd,
                 k = k, rho = rho, lambda = lambda,
                 rho_rule = "direct coefficient, unchanged on full-sample refit"),
      converged = conv, iterations = init_iters + tgd$iterations
    )
  }
  cross_validate_loading_grid(
    full_views, full_prep, fold_objects, rank, grid,
    prepare_context, fit_candidate, label = "SGCA", seed = seed,
    parallel_folds = parallel_folds
  )
}

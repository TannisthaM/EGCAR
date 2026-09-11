# Matrix/block utilities, spectral caches, loading extraction and fit assembly.

egcar_project <- function(z, e, A) {
  k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
  if (z$project_left[[e]]) crossprod(z$Q[[k]], A) %*% z$Q[[l]] else
    crossprod(z$Q[[k]], A %*% z$Q[[l]])
}

egcar_lift <- function(z, e, A) {
  k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
  if (z$lift_left[[e]]) tcrossprod(z$Q[[k]] %*% A, z$Q[[l]]) else
    z$Q[[k]] %*% tcrossprod(A, z$Q[[l]])
}

egcar_prepare_context <- function(prep, centered_views = NULL) {
  spectra <- lapply(seq_len(prep$K), function(k) {
    if (!is.null(centered_views) && nrow(centered_views[[k]]) < prep$p_list[[k]]) {
      # Statistical covariance is still crossprod(X)/n. Thin SVD avoids
      # propagating numerical eigenvalues in the algebraic covariance nullspace.
      X <- centered_views[[k]]
      ss <- svd(X, nu = 0L, nv = min(dim(X)))
      keep <- which(ss$d > 0)
      return(list(vectors = ss$v[, keep, drop = FALSE], values = ss$d[keep]^2 / nrow(X),
                  origin = "thin training-data SVD"))
    }
    ev <- prep$eig[[k]]; keep <- which(ev$values > 0)
    list(vectors = ev$vectors[, keep, drop = FALSE], values = ev$values[keep],
         origin = "cached covariance eigendecomposition")
  })
  Q <- lapply(spectra, `[[`, "vectors")
  d <- lapply(spectra, `[[`, "values")
  rk <- lengths(d)
  ek <- prep$edge_table$k; el <- prep$edge_table$l
  pk <- prep$p_list[ek]; pl <- prep$p_list[el]
  a <- rk[ek]; b <- rk[el]
  # Compare total multiply counts, not just the size of one intermediate.
  z <- list(Q = Q, p_list = as.integer(prep$p_list), edge_k = as.integer(ek),
    edge_l = as.integer(el), keys = prep$edge_table$key, q = prep$q,
    cols_k = prep$edge_cols_k, cols_l = prep$edge_cols_l,
    S = unname(prep$S_kl), full = (a == pk & b == pl),
    basis_origin = vapply(spectra, `[[`, character(1L), "origin"),
    project_left = (a * pk * pl + a * pl * b <= pk * pl * b + a * pk * b),
    lift_left = (pk * a * b + pk * b * pl <= a * b * pl + pk * a * pl))
  z$D <- lapply(seq_along(ek), function(e) outer(d[[ek[[e]]]], d[[el[[e]]]], "*"))
  z$St <- lapply(seq_along(ek), function(e) egcar_project(z, e, z$S[[e]]))
  z$remainder <- lapply(seq_along(ek), function(e) {
    if (z$full[[e]]) matrix(0, 0L, 0L) else z$S[[e]] - egcar_lift(z, e, z$St[[e]])
  })
  z
}

cache_problem_matrices <- function(prep) {
  prep <- cache_problem_matrices_reference(prep)
  # The data-aware prepare_problem wrapper handles wide training views.
  if (EGCAR_BACKEND != "reference" && !any(prep$p_list > prep$n))
    prep$egcar_context <- egcar_prepare_context(prep)
  prep
}

prepare_problem <- function(centered_views) {
  prep <- prepare_problem_reference(centered_views)
  if (EGCAR_BACKEND != "reference" && is.null(prep$egcar_context))
    prep$egcar_context <- egcar_prepare_context(prep, centered_views)
  prep
}

egcar_get_context <- function(prep) {
  if (!is.null(prep$egcar_context)) return(prep$egcar_context)
  if (is.null(prep$edge_cols_k)) prep <- cache_problem_matrices_reference(prep)
  egcar_prepare_context(prep)
}

egcar_initial_state <- function(prep, z, init, group) {
  C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
  # Match the original default initialization of every consensus copy.
  if (!group) return(list(C = C,
    Z = if (!is.null(init$Z)) init$Z else empty_edge_list(prep$edge_table),
    H = if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)))
  G <- if (!is.null(init$G)) init$G else assemble_all_M(C, prep$layout)
  V <- if (!is.null(init$V)) init$V else lapply(G, function(A) matrix(0, nrow(A), ncol(A)))
  E <- seq_along(z$edge_k)
  list(C = C,
    Gk = lapply(E, function(e) G[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
    Gl = lapply(E, function(e) t(G[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])),
    Vk = lapply(E, function(e) V[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
    Vl = lapply(E, function(e) t(V[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])))
}

egcar_loading_factors <- function(prep, selected, covariance_ridge) {
  key <- paste0("egcar-block:", sprintf("%.17g", covariance_ridge), ":", paste(selected, collapse = ","))
  cache <- prep$loading_factor_cache
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE))
    return(get(key, envir = cache, inherits = FALSE))
  idx <- prep$indices %||% make_block_indices(prep$p_list)
  local <- lapply(idx, function(ii) which(ii %in% selected))
  position <- lapply(seq_along(idx), function(k) match(idx[[k]][local[[k]]], selected))
  diag_selected <- unlist(lapply(seq_along(idx), function(k) diag(prep$S_kk[[k]])[local[[k]]]), use.names = FALSE)
  scale_diag <- mean(diag_selected)
  if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
  blocks <- lapply(seq_along(idx), function(k) {
    sk <- local[[k]]
    if (!length(sk)) return(list(half = matrix(0, 0L, 0L), inv_half = matrix(0, 0L, 0L)))
    ev <- eigen(symmetrize(prep$S_kk[[k]][sk, sk, drop = FALSE]), symmetric = TRUE)
    d <- pmax(ev$values + covariance_ridge * scale_diag, 1e-10)
    list(half = tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors),
         inv_half = tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors))
  })
  ans <- list(local = local, position = position, blocks = blocks)
  if (is.environment(cache) && length(cache) < LOADING_FACTOR_CACHE_MAX) assign(key, ans, envir = cache)
  ans
}

egcar_top_eigen <- function(A, rank, positive_tol = 1e-10) {
  n <- nrow(A)
  # Full eigen is usually cheaper for the small default p=45 experiment.
  use_partial <- EGCAR_PARTIAL_EIGEN && n >= EGCAR_PARTIAL_EIGEN_MIN &&
    4L * rank < n && requireNamespace("RSpectra", quietly = TRUE)
  if (use_partial) {
    ee <- tryCatch(suppressWarnings(RSpectra::eigs_sym(A, k = rank, which = "LA",
      opts = list(tol = 1e-12, maxitr = 2000L,
                  ncv = min(n, max(4L * rank + 1L, 30L)),
                  initvec = sin(seq_len(n)) + cos(sqrt(2) * seq_len(n))))),
      error = function(e) NULL)
    if (!is.null(ee) && ee$nconv == rank && length(ee$values) == rank &&
        all(is.finite(ee$values)) && all(is.finite(ee$vectors))) {
      ord <- order(ee$values, decreasing = TRUE)
      ee$values <- ee$values[ord]; ee$vectors <- ee$vectors[, ord, drop = FALSE]
      residual <- A %*% ee$vectors - sweep(ee$vectors, 2L, ee$values, "*")
      error <- sqrt(colSums(residual * residual))
      scale <- max(1, max(rowSums(abs(A))))
      orth_error <- max(abs(crossprod(ee$vectors) - diag(rank)))
      # Resolve numerical rank/positivity near the original threshold using
      # the original dense solver, rather than letting iterative error decide.
      away_from_cutoff <- abs(ee$values[[rank]] - positive_tol) >
        max(1e-8 * scale, 4 * max(error))
      if (max(error) <= 1e-9 * scale && orth_error <= 1e-8 && away_from_cutoff)
        return(list(values = ee$values, vectors = ee$vectors))
    }
  }
  ee <- eigen(A, symmetric = TRUE)
  list(values = ee$values[seq_len(rank)], vectors = ee$vectors[, seq_len(rank), drop = FALSE])
}

egcar_loading_from_operator <- function(prep, C, rank, row_threshold = 1e-4,
    covariance_ridge = 1e-4, require_positive = TRUE, positive_tol = 1e-10,
    keep_full_C = TRUE) {
  if (EGCAR_BACKEND == "reference") return(loading_from_operator(
    prep, C, rank, row_threshold, covariance_ridge, require_positive, positive_tol))
  z <- egcar_get_context(prep)
  norms <- egcar_group_norms(z, C)
  selected <- which(unlist(norms, use.names = FALSE) > row_threshold)
  if (length(selected) < rank) return(list(valid = FALSE, reason = "fewer selected rows than rank"))
  f <- egcar_loading_factors(prep, selected, covariance_ridge)
  Rsel <- matrix(0, length(selected), length(selected))
  for (e in seq_along(z$edge_k)) {
    k <- z$edge_k[[e]]; l <- z$edge_l[[e]]
    sk <- f$local[[k]]; sl <- f$local[[l]]
    if (!length(sk) || !length(sl)) next
    B <- f$blocks[[k]]$half %*% C[[e]][sk, sl, drop = FALSE] %*% f$blocks[[l]]$half
    Rsel[f$position[[k]], f$position[[l]]] <- B
    Rsel[f$position[[l]], f$position[[k]]] <- t(B)
  }
  ee <- egcar_top_eigen(symmetrize(Rsel), rank, positive_tol)
  if (require_positive && ee$values[[rank]] <= positive_tol)
    return(list(valid = FALSE, reason = "fewer than rank positive eigenvalues"))
  U <- ee$vectors
  L <- matrix(0, prep$p, rank)
  for (k in seq_along(z$p_list)) {
    ii <- f$position[[k]]
    if (length(ii)) L[selected[ii], ] <- f$blocks[[k]]$inv_half %*% U[ii, , drop = FALSE]
  }
  list(valid = TRUE, L = L, U = U, selected = selected, eigenvalues = ee$values,
    generalized_eigenvalues = 1 + ee$values,
    C_full = if (keep_full_C) assemble_full_C(C, prep$p_list) else NULL)
}

symmetrize <- function(A) {
  (A + t(A)) / 2
}

frob <- function(A) {
  sqrt(sum(A * A))
}

row_l2 <- function(A) {
  sqrt(rowSums(A * A))
}

soft_threshold <- function(A, threshold) {
  if (threshold <= 0) return(A)
  sign(A) * pmax(abs(A) - threshold, 0)
}

row_group_threshold <- function(A, threshold) {
  if (threshold <= 0) return(A)
  nr <- row_l2(A)
  mult <- pmax(0, 1 - threshold / pmax(nr, .Machine$double.eps))
  A * mult
}

matrix_power_psd <- function(A, power, ridge = 0, eig_floor = 1e-10) {
  ev <- eigen(symmetrize(A), symmetric = TRUE)
  values <- pmax(ev$values + ridge, eig_floor)
  tcrossprod(sweep(ev$vectors, 2L, values^power, "*"), ev$vectors)
}

block_diag <- function(blocks) {
  nr <- sum(vapply(blocks, nrow, integer(1L)))
  nc <- sum(vapply(blocks, ncol, integer(1L)))
  out <- matrix(0, nr, nc)
  r0 <- 1L
  c0 <- 1L
  for (B in blocks) {
    rr <- r0:(r0 + nrow(B) - 1L)
    cc <- c0:(c0 + ncol(B) - 1L)
    out[rr, cc] <- B
    r0 <- max(rr) + 1L
    c0 <- max(cc) + 1L
  }
  out
}

make_block_indices <- function(p_list) {
  starts <- cumsum(c(1L, head(p_list, -1L)))
  lapply(seq_along(p_list), function(k) {
    seq.int(starts[[k]], length.out = p_list[[k]])
  })
}

edge_key <- function(k, l) {
  if (k > l) {
    tmp <- k
    k <- l
    l <- tmp
  }
  paste0(k, "_", l)
}

make_edge_table <- function(p_list) {
  K <- length(p_list)
  ij <- which(upper.tri(matrix(FALSE, K, K)), arr.ind = TRUE)
  data.frame(
    k = ij[, 1L],
    l = ij[, 2L],
    key = apply(ij, 1L, function(z) edge_key(z[[1L]], z[[2L]])),
    p_k = p_list[ij[, 1L]],
    p_l = p_list[ij[, 2L]],
    stringsAsFactors = FALSE
  )
}

make_incidence_layout <- function(p_list) {
  K <- length(p_list)
  edges <- make_edge_table(p_list)
  lapply(seq_len(K), function(k) {
    neighbors <- setdiff(seq_len(K), k)
    widths <- p_list[neighbors]
    ends <- cumsum(widths)
    starts <- c(1L, head(ends, -1L) + 1L)
    cols <- lapply(seq_along(neighbors), function(j) starts[[j]]:ends[[j]])
    names(cols) <- as.character(neighbors)
    list(
      view = k,
      edge_ids = vapply(neighbors, function(l) {
        which(edges$k == min(k, l) & edges$l == max(k, l))
      }, integer(1L)),
      neighbors = neighbors,
      cols = cols,
      nrow = p_list[[k]],
      ncol = sum(widths)
    )
  })
}

empty_edge_list <- function(edge_table, value = 0) {
  out <- lapply(seq_len(nrow(edge_table)), function(e) {
    matrix(value, edge_table$p_k[[e]], edge_table$p_l[[e]])
  })
  names(out) <- edge_table$key
  out
}

assemble_Mk <- function(C, k, layout) {
  lk <- layout[[k]]
  pieces <- lapply(seq_along(lk$neighbors), function(j) {
    l <- lk$neighbors[[j]]
    A <- if (!is.null(lk$edge_ids)) C[[lk$edge_ids[[j]]]] else C[[edge_key(k, l)]]
    if (k < l) A else t(A)
  })
  do.call(cbind, pieces)
}

assemble_all_M <- function(C, layout) {
  lapply(seq_along(layout), function(k) assemble_Mk(C, k, layout))
}

extract_edge_slice <- function(W, k, l, endpoint, layout) {
  stopifnot(k < l, endpoint %in% c(k, l))
  other <- if (endpoint == k) l else k
  cc <- layout[[endpoint]]$cols[[as.character(other)]]
  block <- W[, cc, drop = FALSE]
  if (endpoint == k) block else t(block)
}

assemble_full_C <- function(C, p_list) {
  idx <- make_block_indices(p_list)
  edge_table <- make_edge_table(p_list)
  p <- sum(p_list)
  out <- matrix(0, p, p)
  for (e in seq_len(nrow(edge_table))) {
    k <- edge_table$k[[e]]
    l <- edge_table$l[[e]]
    A <- C[[edge_table$key[[e]]]]
    out[idx[[k]], idx[[l]]] <- A
    out[idx[[l]], idx[[k]]] <- t(A)
  }
  out
}

split_full_C <- function(C_full, p_list) {
  idx <- make_block_indices(p_list)
  edge_table <- make_edge_table(p_list)
  out <- setNames(vector("list", nrow(edge_table)), edge_table$key)
  for (e in seq_len(nrow(edge_table))) {
    k <- edge_table$k[[e]]
    l <- edge_table$l[[e]]
    out[[edge_table$key[[e]]]] <-
      C_full[idx[[k]], idx[[l]], drop = FALSE]
  }
  out
}

center_views <- function(views) {
  means <- lapply(views, colMeans)
  centered <- Map(function(X, m) sweep(X, 2L, m, "-"), views, means)
  list(views = centered, means = means)
}

center_views_at <- function(views, means) {
  Map(function(X, m) sweep(X, 2L, m, "-"), views, means)
}

pairwise_loss <- function(C, S_kk, S_ll, S_kl) {
  0.5 * sum(C * (S_kk %*% C %*% S_ll)) - sum(S_kl * C)
}

operator_objective <- function(prep, C, rho_e, lambda_g) {
  loss <- 0
  for (e in seq_len(nrow(prep$edge_table))) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    key <- prep$edge_table$key[[e]]
    loss <- loss + pairwise_loss(
      C[[key]], prep$S_kk[[k]], prep$S_kk[[l]], prep$S_kl[[key]]
    )
  }
  M <- assemble_all_M(C, prep$layout)
  loss + rho_e * sum(vapply(C, function(A) sum(abs(A)), numeric(1L))) +
    lambda_g * sum(vapply(M, function(A) sum(row_l2(A)), numeric(1L)))
}

loading_metric_factors <- function(prep, selected, covariance_ridge) {
  key <- paste0(sprintf("%.17g", covariance_ridge), ":", paste(selected, collapse = ","))
  cache <- prep$loading_factor_cache
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  Sigma0 <- prep$Sigma0 %||% block_diag(prep$S_kk)
  S <- Sigma0[selected, selected, drop = FALSE]
  scale_diag <- mean(diag(S))
  if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
  ev <- eigen(symmetrize(S), symmetric = TRUE)
  d <- pmax(ev$values + covariance_ridge * scale_diag, 1e-10)
  ans <- list(
    half = tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors),
    inv_half = tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
  )
  if (is.environment(cache) && length(cache) < LOADING_FACTOR_CACHE_MAX) {
    assign(key, ans, envir = cache)
  }
  ans
}

make_validation_covariance <- function(centered_views) {
  n <- nrow(centered_views[[1L]])
  X <- do.call(cbind, centered_views)
  Sigma <- crossprod(X) / n
  # Diagonal blocks are already present in Sigma; do not traverse X again.
  idx <- make_block_indices(vapply(centered_views, ncol, integer(1L)))
  Sigma0 <- block_diag(lapply(idx, function(ii) Sigma[ii, ii, drop = FALSE]))
  list(Sigma = Sigma, Sigma0 = Sigma0)
}

.egcar_validate_init <- function(init, prep, penalty) {
  if (is.null(init)) return(NULL)
  if (inherits(init, "egcar_fit")) {
    if (!identical(init$penalty, penalty) || !identical(init$p_list, prep$p_list))
      stop("The warm start must use the same penalty family and view dimensions.")
    init <- init$solver
  }
  if (!is.list(init)) stop("init must be a fit or an ADMM-state list.")
  edge_fields <- if (penalty == "l11") c("C", "Z", "H") else "C"
  check <- function(x, nr, nc, field) {
    if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(as.integer(nr), as.integer(nc))) ||
        any(!is.finite(x))) stop("Invalid warm-start matrix in ", field, ".")
  }
  for (nm in edge_fields) if (!is.null(init[[nm]])) {
    if (!is.list(init[[nm]]) || length(init[[nm]]) != nrow(prep$edge_table))
      stop("Invalid warm-start list: ", nm)
    if (!is.null(names(init[[nm]]))) {
      if (!setequal(names(init[[nm]]), prep$edge_table$key)) stop("Warm-start edge names do not match.")
      init[[nm]] <- init[[nm]][prep$edge_table$key]
    }
    for (j in seq_len(nrow(prep$edge_table)))
      check(init[[nm]][[j]], prep$edge_table$p_k[[j]], prep$edge_table$p_l[[j]], nm)
  }
  if (penalty == "l21") for (nm in c("G", "V")) if (!is.null(init[[nm]])) {
    if (!is.list(init[[nm]]) || length(init[[nm]]) != prep$K) stop("Invalid warm-start list: ", nm)
    for (k in seq_len(prep$K)) check(init[[nm]][[k]], prep$p_list[[k]], prep$p - prep$p_list[[k]], nm)
  }
  init
}

.egcar_loading <- function(prep, solver, rank, control, e, keep_full_C = control$keep_full_C) {
  out <- e$egcar_loading_from_operator(prep, solver$C_hat, rank,
    row_threshold = control$row_threshold, covariance_ridge = control$covariance_ridge,
    require_positive = TRUE, keep_full_C = keep_full_C)
  if (!keep_full_C) out$C_full <- NULL
  out
}

.egcar_fit_object <- function(raw, loading, prepared, rank, penalty, lambda,
                              control, fit_time, loading_time, call = NULL) {
  ok <- isTRUE(loading$valid)
  L <- if (ok) loading$L else NULL
  if (!is.null(L)) {
    rownames(L) <- unlist(Map(function(v, nm) paste(v, nm, sep = ":"),
                               prepared$view_names, prepared$feature_names), use.names = FALSE)
    colnames(L) <- paste0("component", seq_len(rank))
  }
  status <- if (!ok) "invalid_loading" else if (!isTRUE(raw$converged)) "not_converged" else "ok"
  structure(list(L = L, C = raw$C_hat, C_full = loading$C_full,
    loadings = if (ok) setNames(lapply(make_block_indices(prepared$p_list),
       function(ii) L[ii, , drop = FALSE]), prepared$view_names) else NULL,
    loading = loading, solver = raw, rank = rank, penalty = penalty, lambda = lambda,
    rho_e = if (penalty == "l11") lambda else 0,
    lambda_g = if (penalty == "l21") lambda else 0,
    selected = if (ok) loading$selected else integer(),
    means = prepared$means, view_names = prepared$view_names, feature_names = prepared$feature_names,
    p_list = prepared$p_list, p = prepared$p, n = prepared$n,
    converged = raw$converged, iterations = raw$iterations,
    status = status, error = if (!ok) loading$reason else NULL,
    control = control, fit_time = fit_time, loading_time = loading_time,
    preparation_time = prepared$preparation_time, call = call), class = "egcar_fit")
}

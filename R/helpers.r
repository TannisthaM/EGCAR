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
         origin = ev$origin %||% "cached covariance eigendecomposition")
  })
  Q <- lapply(spectra, `[[`, "vectors")
  d <- lapply(spectra, `[[`, "values")
  rk <- lengths(d)
  ek <- prep$edge_table$k; el <- prep$edge_table$l
  # pk/pl/a/b feed the triple products below (a*pk*pl etc.). p_list and rk
  # (= lengths(d)) are native R integers, and at p_total in the thousands
  # (e.g. p_per_block=5000, where a population covariance is generically
  # full rank so a/b ~= 5000 too) a*pk*pl alone reaches ~1.25e11 -- well
  # past the 32-bit integer ceiling of ~2.15e9. Plain int*int*int silently
  # overflows to NA (with a warning) rather than erroring, so this must be
  # forced into double arithmetic; doubles represent integers exactly up to
  # 2^53, far beyond anything these dimensions reach.
  pk <- as.double(prep$p_list[ek]); pl <- as.double(prep$p_list[el])
  a <- as.double(rk[ek]); b <- as.double(rk[el])
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
  if (EGCAR_BACKEND == "reference") return(cache_problem_matrices_reference(prep))
  prep$indices <- make_block_indices(prep$p_list)
  prep$edge_cols_k <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
    k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
    prep$layout[[k]]$cols[[as.character(l)]]
  })
  prep$edge_cols_l <- lapply(seq_len(nrow(prep$edge_table)), function(e) {
    k <- prep$edge_table$k[[e]]; l <- prep$edge_table$l[[e]]
    prep$layout[[l]]$cols[[as.character(k)]]
  })
  prep$loading_factor_cache <- new.env(parent = emptyenv())
  prep$egcar_context <- egcar_prepare_context(prep)
  prep
}

prepare_problem <- function(centered_views) {
  if (EGCAR_BACKEND == "reference") return(prepare_problem_reference(centered_views))
  n <- nrow(centered_views[[1L]])
  p_list <- vapply(centered_views, ncol, integer(1L))
  if (any(vapply(centered_views, nrow, integer(1L)) != n))
    stop("All views must have the same number of rows.")
  edges <- make_edge_table(p_list)
  S_kk <- lapply(centered_views, function(X) crossprod(X) / n)
  S_kl <- setNames(lapply(seq_len(nrow(edges)), function(e)
    crossprod(centered_views[[edges$k[[e]]]], centered_views[[edges$l[[e]]]]) / n), edges$key)
  eig <- lapply(seq_along(p_list), function(k) {
    if (n < p_list[[k]]) {
      # Use the same thin basis as before, without first computing and storing
      # an unused p_k by p_k covariance eigendecomposition. No new cutoff.
      ss <- svd(centered_views[[k]], nu = 0L, nv = min(n, p_list[[k]]))
      keep <- which(ss$d > 0)
      list(vectors = ss$v[, keep, drop = FALSE], values = ss$d[keep]^2 / n,
           origin = "thin training-data SVD")
    } else {
      ee <- eigen(symmetrize(S_kk[[k]]), symmetric = TRUE)
      ee$values <- pmax(ee$values, 0)
      ee
    }
  })
  cache_problem_matrices(list(n = n, K = length(p_list), p_list = p_list,
    p = sum(p_list), edge_table = edges, layout = make_incidence_layout(p_list),
    S_kk = S_kk, S_kl = S_kl, eig = eig,
    q = sum(as.double(edges$p_k) * edges$p_l)))
}

egcar_get_context <- function(prep) {
  if (!is.null(prep$egcar_context)) return(prep$egcar_context)
  if (is.null(prep$edge_cols_k)) return(cache_problem_matrices(prep)$egcar_context)
  egcar_prepare_context(prep)
}

egcar_initial_state <- function(prep, z, init, group) {
  C <- if (!is.null(init$C)) init$C else empty_edge_list(prep$edge_table)
  # Match the original default initialization of every consensus copy.
  if (!group) return(list(C = C,
    Z = if (!is.null(init$Z)) init$Z else empty_edge_list(prep$edge_table),
    H = if (!is.null(init$H)) init$H else empty_edge_list(prep$edge_table)))
  if (all(c("Hk", "Hl", "a") %in% names(init)))
    return(c(list(C = C), init[c("Hk", "Hl", "a")]))
  endpoint_fields <- c("Gk", "Gl", "Vk", "Vl")
  if (all(endpoint_fields %in% names(init)))
    return(c(list(C = C), init[endpoint_fields]))
  if (is.null(init$G) && is.null(init$V)) {
    # Cold starts need no assembled p_k by (p-p_k) consensus matrices.
    return(list(C = C, Hk = unname(C), Hl = unname(C),
      a = lapply(z$p_list, function(pk) rep.int(1, pk))))
  }
  G <- if (!is.null(init$G)) init$G else assemble_all_M(C, prep$layout)
  V <- if (!is.null(init$V)) init$V else lapply(G, function(A) matrix(0, nrow(A), ncol(A)))
  E <- seq_along(z$edge_k)
  list(C = C,
    Gk = lapply(E, function(e) G[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
    Gl = lapply(E, function(e) t(G[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])),
    Vk = lapply(E, function(e) V[[z$edge_k[[e]]]][, z$cols_k[[e]], drop = FALSE]),
    Vl = lapply(E, function(e) t(V[[z$edge_l[[e]]]][, z$cols_l[[e]], drop = FALSE])))
}

.egcar_cache_key <- function(prefix, covariance_ridge, selected) {
  # R caps every variable/symbol name (including exists()/assign()/get()
  # targets on an environment) at 10000 bytes. The original key here was
  # built from paste(selected, collapse=","), a full comma-joined row-index
  # list -- fine for small p, but at p_total in the low thousands (already
  # observed at p_per_block=1000) enough rows get selected that the joined
  # string exceeds that cap, crashing with "variable names are limited to
  # 10000 bytes". This builds a short, fixed-width surrogate key instead.
  # It is NOT collision-proof by itself (it doesn't need to be): the actual
  # `selected` vector and `covariance_ridge` are stored alongside the cached
  # value and checked with identical() at the call site before a hit is
  # trusted, so a collision only costs a harmless cache miss/recompute, never
  # a wrong cached result being returned silently.
  n <- length(selected)
  fold <- if (n) sum((as.numeric(selected) %% 999983) * seq_len(n)) %% 2147483647 else 0
  sprintf("%s%.17g:%d:%.0f", prefix, covariance_ridge, n, fold)
}

egcar_cache_loading_factors <- function(cache, key, selected, ridge, value) {
  if (!is.environment(cache) || length(cache) >= LOADING_FACTOR_CACHE_MAX ||
      LOADING_FACTOR_CACHE_MAX_BYTES <= 0) return(invisible(NULL))
  entry <- list(sel = selected, ridge = ridge, value = value)
  needed <- as.numeric(utils::object.size(entry))
  if (needed > LOADING_FACTOR_CACHE_MAX_BYTES) return(invisible(NULL))
  keys <- setdiff(ls(cache, all.names = TRUE), key)
  used <- sum(vapply(keys, function(k)
    as.numeric(utils::object.size(get(k, envir = cache, inherits = FALSE))), numeric(1L)))
  if (used + needed <= LOADING_FACTOR_CACHE_MAX_BYTES)
    assign(key, entry, envir = cache)
  invisible(NULL)
}

egcar_loading_factors <- function(prep, selected, covariance_ridge) {
  # New key separates compact factors from pre-0.2.12 dense cache entries.
  key <- .egcar_cache_key("egcar-operator:", covariance_ridge, selected)
  cache <- prep$loading_factor_cache
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
    hit <- get(key, envir = cache, inherits = FALSE)
    if (identical(hit$sel, selected) && identical(hit$ridge, covariance_ridge)) return(hit$value)
  }
  idx <- prep$indices %||% make_block_indices(prep$p_list)
  local <- lapply(idx, function(ii) which(ii %in% selected))
  position <- lapply(seq_along(idx), function(k) match(idx[[k]][local[[k]]], selected))
  diag_selected <- unlist(lapply(seq_along(idx), function(k) diag(prep$S_kk[[k]])[local[[k]]]), use.names = FALSE)
  scale_diag <- mean(diag_selected)
  if (!is.finite(scale_diag) || scale_diag <= 0) scale_diag <- 1
  ridge <- covariance_ridge * scale_diag
  blocks <- lapply(seq_along(idx), function(k) {
    sk <- local[[k]]; old <- prep$eig[[k]]
    if (!length(sk)) return(list(vectors = matrix(0, 0L, 0L),
      half = numeric(), inv_half = numeric(), base_half = 0, base_inv_half = 0))
    if (length(sk) == prep$p_list[[k]]) {
      ev <- old
    } else if (ncol(old$vectors) < length(sk)) {
      # S[sk,sk] = F F'. SVD of this thin F keeps the complete covariance
      # range without forming a selected covariance or a dense square root.
      F <- sweep(old$vectors[sk, , drop = FALSE], 2L, sqrt(pmax(old$values, 0)), "*")
      if (ncol(F)) {
        sv <- svd(F, nu = min(dim(F)), nv = 0L)
        ev <- list(vectors = sv$u, values = sv$d^2)
      } else ev <- list(vectors = matrix(0, length(sk), 0L), values = numeric())
    } else {
      ev <- eigen(symmetrize(prep$S_kk[[k]][sk, sk, drop = FALSE]), symmetric = TRUE)
    }
    # Preserve the established scale-relative ridge and 1e-10 metric floor.
    d <- pmax(ev$values + ridge, 1e-10)
    base <- if (ncol(ev$vectors) < length(sk)) max(ridge, 1e-10) else 0
    bh <- sqrt(base); bi <- if (base > 0) 1 / bh else 0
    list(vectors = ev$vectors, half = sqrt(d) - bh,
      inv_half = 1 / sqrt(d) - bi, base_half = bh, base_inv_half = bi)
  })
  ans <- list(local = local, position = position, blocks = blocks)
  egcar_cache_loading_factors(cache, key, selected, covariance_ridge, ans)
  ans
}

egcar_apply_loading_factor <- function(block, x, inverse = FALSE) {
  base <- if (inverse) block$base_inv_half else block$base_half
  weights <- if (inverse) block$inv_half else block$half
  if (!length(weights)) return(base * x)
  ans <- base * x + block$vectors %*% (weights * crossprod(block$vectors, x))
  if (is.null(dim(x))) as.numeric(ans) else ans
}

egcar_apply_loading_metric <- function(factors, v, inverse = FALSE) {
  ans <- numeric(length(v))
  for (k in seq_along(factors$blocks)) {
    ii <- factors$position[[k]]
    if (length(ii)) ans[ii] <- egcar_apply_loading_factor(factors$blocks[[k]], v[ii], inverse)
  }
  ans
}

egcar_operator_multiply <- function(v, args) {
  if (args$native) return(egcar_native_block_product(args$C, args$edge_k,
    args$edge_l, args$factors$local, args$factors$position, v))
  # Independent R/BLAS path, without repeatedly copying selected edge matrices.
  f <- args$factors
  x <- lapply(seq_along(args$p_list), function(k) {
    out <- numeric(args$p_list[[k]])
    out[f$local[[k]]] <- v[f$position[[k]]]; out
  })
  out <- numeric(length(v))
  for (e in seq_along(args$edge_k)) {
    k <- args$edge_k[[e]]; l <- args$edge_l[[e]]
    if (!length(f$local[[k]]) || !length(f$local[[l]])) next
    out[f$position[[k]]] <- out[f$position[[k]]] +
      as.numeric(args$C[[e]] %*% x[[l]])[f$local[[k]]]
    out[f$position[[l]]] <- out[f$position[[l]]] +
      as.numeric(crossprod(args$C[[e]], x[[k]]))[f$local[[l]]]
  }
  out
}

egcar_loading_operator <- function(v, args) {
  x <- egcar_apply_loading_metric(args$factors, v)
  egcar_apply_loading_metric(args$factors, egcar_operator_multiply(x, args))
}

egcar_top_eigen <- function(A, rank, positive_tol = 1e-10, n = NULL, args = NULL) {
  if (is.null(n)) n <- nrow(A)
  if (length(n) != 1L || rank < 1L || rank > n)
    stop("The requested eigen rank must be between one and the selected dimension.")
  multiply <- if (is.function(A)) function(v) A(v, args) else function(v) as.numeric(A %*% v)
  use_partial <- EGCAR_PARTIAL_EIGEN && n >= EGCAR_PARTIAL_EIGEN_MIN && rank < n
  if (use_partial && n == 2L) {
    # RSpectra requires n >= 3. The leading pair of a 2x2 symmetric matrix
    # has a closed form, so this case needs neither RSpectra nor a full EVD.
    left <- multiply(c(1, 0)); right <- multiply(c(0, 1))
    delta <- left[[1L]] / 2 - right[[2L]] / 2
    off <- left[[2L]] / 2 + right[[1L]] / 2
    scale <- max(abs(delta), abs(off))
    radius <- if (scale == 0) 0 else scale * sqrt((delta / scale)^2 + (off / scale)^2)
    angle <- atan2(off, delta) / 2
    u <- matrix(c(cos(angle), sin(angle)), 2L, 1L)
    value <- left[[1L]] / 2 + right[[2L]] / 2 + radius
    return(list(values = value, vectors = u, method = "analytic 2x2", nconv = 1L,
      residuals = sqrt(sum((multiply(as.numeric(u)) - value * u)^2)), nops = 3L))
  }
  if (use_partial) {
    if (!requireNamespace("RSpectra", quietly = TRUE))
      stop("Matrix-free loading extraction requires RSpectra. Install it and retry.")
    # A larger Krylov space is the only automatic retry. Never silently allocate
    # a dense selected operator or compute a full EVD after nonconvergence.
    last_error <- "no converged eigenpairs"
    for (attempt in 1:2) {
      ee <- tryCatch(RSpectra::eigs_sym(A, k = rank, n = n, args = args, which = "LA",
        opts = list(tol = 1e-12, maxitr = 2000L * attempt,
          ncv = as.integer(min(n, max((4L * attempt) * rank + 1L, 30L * attempt))),
          initvec = sin(seq_len(n)) + cos(sqrt(2) * seq_len(n)))),
        error = function(e) { last_error <<- conditionMessage(e); NULL })
      if (is.null(ee) || ee$nconv != rank || length(ee$values) != rank ||
          any(!is.finite(ee$values)) || any(!is.finite(ee$vectors))) next
      ord <- order(ee$values, decreasing = TRUE)
      values <- ee$values[ord]; vectors <- ee$vectors[, ord, drop = FALSE]
      errors <- vapply(seq_len(rank), function(j)
        sqrt(sum((multiply(vectors[, j]) - values[[j]] * vectors[, j])^2)), numeric(1L))
      scale <- max(1, abs(values))
      orth_error <- max(abs(crossprod(vectors) - diag(rank)))
      if (max(errors) <= 1e-9 * scale && orth_error <= 1e-8)
        return(list(values = values, vectors = vectors, method = "matrix-free partial",
          residuals = errors, nconv = ee$nconv, nops = ee$nops))
      last_error <- "eigenpair residual or orthogonality check failed"
    }
    stop("Partial loading eigensolve failed: ", last_error,
      ". No dense fallback was attempted. For a small diagnostic problem, explicitly set partial_eigen = FALSE.")
  }
  # Explicit diagnostic mode, a user-selected minimum dimension, or rank == n.
  if (is.function(A)) A <- vapply(seq_len(n), function(j) {
    v <- numeric(n); v[[j]] <- 1; multiply(v)
  }, numeric(n))
  ee <- eigen(symmetrize(A), symmetric = TRUE)
  list(values = ee$values[seq_len(rank)], vectors = ee$vectors[, seq_len(rank), drop = FALSE],
    method = "dense", nconv = rank)
}

egcar_loading_from_operator <- function(prep, C, rank, row_threshold = 1e-4,
    covariance_ridge = 1e-4, require_positive = TRUE, positive_tol = 1e-10,
    keep_full_C = FALSE) {
  if (EGCAR_BACKEND == "reference") return(loading_from_operator(
    prep, C, rank, row_threshold, covariance_ridge, require_positive, positive_tol))
  z <- egcar_get_context(prep)
  norms <- egcar_group_norms(z, C)
  selected <- which(unlist(norms, use.names = FALSE) > row_threshold)
  if (length(selected) < rank) return(list(valid = FALSE, reason = "fewer selected rows than rank"))
  f <- egcar_loading_factors(prep, selected, covariance_ridge)
  args <- list(C = unname(C), edge_k = z$edge_k, edge_l = z$edge_l,
    p_list = z$p_list, factors = f, native = EGCAR_BACKEND == "cpp")
  ee <- egcar_top_eigen(egcar_loading_operator, rank, positive_tol,
    n = length(selected), args = args)
  if (require_positive && ee$values[[rank]] <= positive_tol)
    return(list(valid = FALSE, reason = "fewer than rank positive eigenvalues"))
  if (require_positive && !is.null(ee$residuals) &&
      ee$values[[rank]] - ee$residuals[[rank]] <= positive_tol)
    return(list(valid = FALSE, reason = "positive eigenvalue threshold unresolved at eigensolver precision"))
  U <- ee$vectors
  L <- matrix(0, prep$p, rank)
  for (k in seq_along(z$p_list)) {
    ii <- f$position[[k]]
    if (length(ii)) L[selected[ii], ] <- egcar_apply_loading_factor(f$blocks[[k]],
      U[ii, , drop = FALSE], inverse = TRUE)
  }
  list(valid = TRUE, L = L, U = U, selected = selected, eigenvalues = ee$values,
    generalized_eigenvalues = 1 + ee$values, eigen_method = ee$method,
    eigen_residuals = ee$residuals, eigen_nconv = ee$nconv, eigen_nops = ee$nops,
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
  key <- .egcar_cache_key("", covariance_ridge, selected)
  cache <- prep$loading_factor_cache
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
    hit <- get(key, envir = cache, inherits = FALSE)
    if (identical(hit$sel, selected) && identical(hit$ridge, covariance_ridge)) return(hit$value)
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
  egcar_cache_loading_factors(cache, key, selected, covariance_ridge, ans)
  ans
}

make_validation_covariance <- function(centered_views) {
  n <- nrow(centered_views[[1L]])
  p_list <- vapply(centered_views, ncol, integer(1L))
  # For a wide validation split, retain X rather than two p by p matrices.
  # validation_score contracts through the n by rank component scores.
  if (n < sum(p_list)) return(list(views = centered_views, n = n,
    indices = make_block_indices(p_list)))
  X <- do.call(cbind, centered_views)
  Sigma <- crossprod(X) / n
  # Diagonal blocks are already present in Sigma; do not traverse X again.
  idx <- make_block_indices(vapply(centered_views, ncol, integer(1L)))
  Sigma0 <- block_diag(lapply(idx, function(ii) Sigma[ii, ii, drop = FALSE]))
  # Keep the small mean vector so SGCA can form test-centered covariance
  # scores while EGCAR retains its existing training-centered score.
  list(Sigma = Sigma, Sigma0 = Sigma0, mean = colMeans(X), n = n)
}

.egcar_validate_init <- function(init, prep, penalty) {
  if (is.null(init)) return(NULL)
  if (inherits(init, "egcar_fit")) {
    if (!identical(init$penalty, penalty) || !identical(init$p_list, prep$p_list))
      stop("The warm start must use the same penalty family and view dimensions.")
    init <- init$solver
  }
  if (!is.list(init)) stop("init must be a fit or an ADMM-state list.")
  endpoints <- c("Gk", "Gl", "Vk", "Vl")
  compressed <- c("Hk", "Hl", "a")
  has_compressed <- any(compressed %in% names(init))
  has_endpoints <- any(endpoints %in% names(init))
  if (penalty == "l21" && has_compressed) {
    if (!all(compressed %in% names(init)))
      stop("A compressed group warm start must contain Hk, Hl and a.")
    if (has_endpoints || any(c("G", "V") %in% names(init)))
      stop("Supply either compressed or legacy group state, not both.")
    if (!is.list(init$a) || length(init$a) != prep$K)
      stop("Invalid group row-multiplier list a.")
    for (k in seq_len(prep$K))
      if (!is.numeric(init$a[[k]]) || !is.null(dim(init$a[[k]])) ||
          length(init$a[[k]]) != prep$p_list[[k]] ||
          any(!is.finite(init$a[[k]])) || any(init$a[[k]] < 0 | init$a[[k]] > 1))
        stop("Invalid group row multipliers in a.")
  }
  if (penalty == "l21" && has_endpoints && !all(endpoints %in% names(init)))
    stop("A compact group warm start must contain Gk, Gl, Vk and Vl.")
  edge_fields <- if (penalty == "l11") c("C", "Z", "H") else
    c("C", if (has_compressed) c("Hk", "Hl") else if (has_endpoints) endpoints)
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

# Keep compressed endpoint state between group CV candidates and path fits.
egcar_warm_state <- function(raw, group) {
  fields <- if (!group) c("C", "Z", "H") else
    if (all(c("Hk", "Hl", "a") %in% names(raw))) c("C", "Hk", "Hl", "a") else
    if (all(c("Gk", "Gl", "Vk", "Vl") %in% names(raw)))
      c("C", "Gk", "Gl", "Vk", "Vl") else c("C", "G", "V")
  raw[fields]
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
    control = control, solver_time = fit_time,
    fit_time = fit_time + loading_time, loading_time = loading_time,
    total_time = fit_time + loading_time, timing_version = "0.2.15",
    preparation_time = prepared$preparation_time, call = call), class = "egcar_fit")
}

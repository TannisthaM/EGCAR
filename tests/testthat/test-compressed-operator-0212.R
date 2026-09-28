pv <- function(name) getFromNamespace(name, "egcar")

test_that("compressed ADMM matches independent dense updates, including adaptive mu and legacy starts", {
  set.seed(21201)
  for (wide in c(FALSE, TRUE)) {
    n <- if (wide) 7L else 25L
    dims <- if (wide) c(11L, 1L, 9L) else c(4L, 3L, 5L)
    X <- lapply(dims, function(p) matrix(rnorm(n * p), n, p))
    if (wide) X[[1L]][, 1L] <- 0
    prep <- egcar_prepare(X)$prep
    z <- prep$egcar_context
    legacy <- list(C = setNames(lapply(z$S, function(A) A * 0), z$keys))
    for (nm in c("Gk", "Gl", "Vk", "Vl"))
      legacy[[nm]] <- lapply(z$S, function(A) matrix(rnorm(length(A), sd = .02), nrow(A), ncol(A)))
    for (lambda in c(0, .03, 2)) for (adaptive in c(FALSE, TRUE)) {
      args <- list(prep = prep, lambda_g = lambda, max_iter = 31L,
        abs_tol = 1e-10, rel_tol = 1e-10, adaptive_mu = adaptive, balance_ratio = 1.1,
        adapt_every = 2L, check_every = 2L, keep_history = TRUE)
      ref <- pv(".egcar_engine")(egcar_control(backend = "reference"))
      for (init in list(NULL, legacy)) {
        args$init <- init
        before <- serialize(init, NULL)
        expected <- do.call(ref$fit_l21_admm, args)
        for (backend in c("cpp", "R")) {
          e <- pv(".egcar_engine")(egcar_control(backend = backend))
          actual <- do.call(e$fit_l21_admm, c(args, list(compact_state = TRUE)))
          expect_true(all(c("C", "Hk", "Hl", "a") %in% names(actual)))
          expect_false(any(c("G", "V", "Gk", "Gl", "Vk", "Vl") %in% names(actual)))
          old <- e$egcar_expand_group_state(z, actual)
          expect_equal(lapply(actual$C_hat, unname), lapply(expected$C_hat, unname), tolerance = 2e-8)
          expect_equal(lapply(old$C, unname), lapply(expected$C, unname), tolerance = 2e-8)
          expect_equal(e$egcar_view_copies(z, old$Gk, old$Gl), expected$G, tolerance = 2e-8)
          expect_equal(e$egcar_view_copies(z, old$Vk, old$Vl), expected$V, tolerance = 2e-8)
          expect_equal(actual$history, expected$history, tolerance = 2e-8)
          expect_equal(actual$mu_g, expected$mu_g)
        }
        expect_identical(serialize(init, NULL), before)
      }
    }
  }
})

test_that("public compressed state supports cross-backend warm starts and saves two arrays", {
  sim <- egcar_simulate(n = 18, p_list = c(17, 12, 9), active_per_view = 3, seed = 21202)
  prep <- egcar_prepare(sim$views)
  ctl <- egcar_control(max_iter = 35L, abs_tol = 0, rel_tol = 0)
  fit <- egcar_fit(prep, 1L, "l21", .01, ctl)
  expect_true(all(c("Hk", "Hl", "a") %in% names(fit$solver)))
  q <- prep$prep$q
  s <- fit$solver
  expect_equal(sum(lengths(s$C)) + sum(lengths(s$Hk)) + sum(lengths(s$Hl)), 3 * q)
  expect_equal(sum(lengths(s$a)), prep$p)
  before <- serialize(s, NULL)
  ctl$compact_state <- FALSE
  old <- egcar_fit(prep, 1L, "l21", .01, ctl)
  expect_true(all(c("G", "V") %in% names(old$solver)))
  for (backend in c("cpp", "R", "reference")) {
    ctl$backend <- backend
    a <- egcar_fit(prep, 1L, "l21", .005, ctl, init = fit)
    b <- egcar_fit(prep, 1L, "l21", .005, ctl, init = old)
    expect_equal(a$C, b$C, tolerance = 1e-8)
    expect_equal(a$solver$G, b$solver$G, tolerance = 1e-8)
    expect_equal(a$solver$V, b$solver$V, tolerance = 1e-8)
  }
  expect_identical(serialize(fit$solver, NULL), before)
  bad <- fit$solver; bad$a[[1L]][[1L]] <- 2
  expect_error(egcar_fit(prep, 1L, "l21", .01, init = bad), "row multipliers")
  bad <- fit$solver; bad$Hl <- NULL
  expect_error(egcar_fit(prep, 1L, "l21", .01, init = bad), "Hk, Hl and a")
})

test_that("compact metric actions and native selected products match dense algebra", {
  set.seed(21203)
  for (n in c(6L, 40L)) {
    prep <- egcar_prepare(list(matrix(rnorm(n * 14), n, 14),
      matrix(rnorm(n * 11), n, 11), matrix(rnorm(n * 4), n, 4)))$prep
    e <- pv(".egcar_engine")(egcar_control())
    C <- lapply(prep$S_kl, function(A) matrix(rnorm(length(A)), nrow(A), ncol(A)))
    full <- e$assemble_full_C(C, prep$p_list)
    for (selected in list(seq_len(prep$p), c(1:10, 16:24), c(2L, 15L))) {
      f <- e$egcar_loading_factors(prep, selected, 1e-4)
      Sigma0 <- e$block_diag(prep$S_kk)[selected, selected, drop = FALSE]
      d <- eigen(e$symmetrize(Sigma0), symmetric = TRUE)
      values <- pmax(d$values + 1e-4 * mean(diag(Sigma0)), 1e-10)
      root <- tcrossprod(sweep(d$vectors, 2L, sqrt(values), "*"), d$vectors)
      inverse <- tcrossprod(sweep(d$vectors, 2L, 1 / sqrt(values), "*"), d$vectors)
      v <- rnorm(length(selected))
      expect_equal(e$egcar_apply_loading_metric(f, v), as.numeric(root %*% v), tolerance = 1e-8)
      expect_equal(e$egcar_apply_loading_metric(f, v, TRUE), as.numeric(inverse %*% v), tolerance = 1e-7)
      args <- list(C = C, edge_k = prep$edge_table$k, edge_l = prep$edge_table$l,
        p_list = prep$p_list, factors = f, native = TRUE)
      before <- serialize(C, NULL)
      expected <- root %*% full[selected, selected, drop = FALSE] %*% root %*% v
      expect_equal(e$egcar_loading_operator(v, args), as.numeric(expected), tolerance = 1e-8)
      args$native <- FALSE
      expect_equal(e$egcar_loading_operator(v, args), as.numeric(expected), tolerance = 1e-8)
      expect_identical(serialize(C, NULL), before)
      expect_true(all(vapply(f$blocks, function(b) is.null(dim(b$half)) && is.null(dim(b$inv_half)), logical(1))))
    }
  }
})

test_that("partial eigenpairs select largest algebraic values and never silently fall back", {
  e <- pv(".egcar_engine")(egcar_control())
  A <- diag(c(5, 4, 2, 1, -100, -20))
  op <- function(v, args) as.numeric(args %*% v)
  a <- e$egcar_top_eigen(op, 2L, n = 6L, args = A)
  expect_equal(a$values, c(5, 4), tolerance = 1e-10)
  expect_identical(a$method, "matrix-free partial")
  expect_lte(max(a$residuals), 1e-9)
  expect_equal(a$vectors %*% t(a$vectors), diag(c(1, 1, 0, 0, 0, 0)), tolerance = 1e-9)
  expect_error(e$egcar_top_eigen(function(v, args) stop("synthetic failure"), 1L, n = 6L),
    "No dense fallback was attempted")
  all <- e$egcar_top_eigen(op, 6L, n = 6L, args = A)
  expect_equal(all$values, sort(diag(A), decreasing = TRUE))
  small <- e$egcar_top_eigen(op, 1L, n = 2L, args = matrix(c(0, -3, -3, 0), 2L))
  expect_equal(small$values, 3)
  expect_identical(small$method, "analytic 2x2")
  expect_lte(small$residuals, 1e-12)
})

test_that("both penalties use the matrix-free loading path by default", {
  sim <- egcar_simulate(n = 25L, p_list = c(20L, 15L, 12L), rank = 2L,
    active_per_view = 5L, signal = .8, seed = 21204)
  prep <- egcar_prepare(sim$views)
  for (penalty in c("l11", "l21")) for (rank in c(1L, 2L, 5L)) {
    ctl <- egcar_control(max_iter = 60L, abs_tol = 0, rel_tol = 0)
    a <- egcar_fit(prep, rank, penalty, .01, ctl)
    ctl$partial_eigen <- FALSE
    b <- egcar_fit(prep, rank, penalty, .01, ctl)
    expect_identical(a$loading$eigen_method, "matrix-free partial")
    expect_null(a$C_full)
    expect_equal(a$loading$eigenvalues, b$loading$eigenvalues, tolerance = 1e-8)
    expect_equal(tcrossprod(a$loading$U), tcrossprod(b$loading$U), tolerance = 1e-7)
    expect_equal(a$C, b$C, tolerance = 0)
  }
  e <- pv(".egcar_engine")(egcar_control())
  e$assemble_full_C <- function(...) stop("full C must not be assembled")
  out <- e$egcar_loading_from_operator(prep$prep, a$C, 1L, keep_full_C = FALSE)
  expect_true(out$valid)
  expect_identical(out$eigen_method, "matrix-free partial")
})

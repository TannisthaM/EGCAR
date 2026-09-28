# Numerical checks for changes in 0.2.11. These exercise independent dense
# formulas and the retained reference backend, including rank-deficient data.
private <- function(name) getFromNamespace(name, "egcar")

test_that("compact preparation preserves covariances and complete nullspace updates", {
  set.seed(211)
  cases <- list(list(n = 24L, p = c(3L, 5L, 4L)),
                list(n = 9L, p = c(15L, 4L, 12L)),
                list(n = 8L, p = c(1L, 12L)))
  for (case in cases) {
    X <- lapply(case$p, function(p) matrix(rnorm(case$n * p), case$n, p))
    if (case$p[[1L]] > 1L) X[[1L]][, 1L] <- 0
    centered <- private("center_views")(private(".egcar_views")(X))$views
    prep <- egcar_prepare(X)
    ref <- private("prepare_problem_reference")(centered)
    expect_null(prep$prep$Sigma0)
    expect_null(prep$prep$eig_products)
    expect_equal(prep$prep$S_kk, ref$S_kk, tolerance = 1e-12)
    expect_equal(prep$prep$S_kl, ref$S_kl, tolerance = 1e-12)
    for (k in seq_along(X)) {
      ev <- prep$prep$eig[[k]]
      reconstructed <- tcrossprod(sweep(ev$vectors, 2L, ev$values, "*"), ev$vectors)
      expect_equal(unname(reconstructed), unname(ref$S_kk[[k]]), tolerance = 1e-10)
    }
    before <- serialize(prep$prep$egcar_context, NULL)
    for (penalty in c("l11", "l21")) {
      ctl <- egcar_control(max_iter = 45L, abs_tol = 0, rel_tol = 0,
        check_every = 5L, partial_eigen = FALSE, keep_history = penalty == "l21")
      a <- egcar_fit(prep, 1L, penalty, .02, ctl)
      ctl$backend <- "reference"
      b <- egcar_fit(prep, 1L, penalty, .02, ctl)
      expect_equal(unlist(a$C), unlist(b$C), tolerance = 1e-7)
      expect_equal(a$solver$primal_residual, b$solver$primal_residual, tolerance = 1e-7)
      expect_equal(a$solver$dual_residual, b$solver$dual_residual, tolerance = 1e-7)
      expect_identical(a$iterations, b$iterations)
    }
    expect_identical(serialize(prep$prep$egcar_context, NULL), before)
  }
})

test_that("data-space validation equals the original dense score", {
  set.seed(212)
  for (n in c(1L, 4L, 40L)) for (rank in c(1L, 2L)) {
    X <- list(matrix(rnorm(n * 3L) + 10, n, 3L),
              matrix(rnorm(n * 4L) - 5, n, 4L))
    joined <- do.call(cbind, X)
    dense <- list(Sigma = crossprod(joined) / n,
      Sigma0 = private("block_diag")(lapply(X, function(z) crossprod(z) / n)))
    compact <- private("make_validation_covariance")(X)
    L <- matrix(rnorm(7L * rank), 7L, rank)
    score <- private("validation_score")
    expect_equal(score(L, compact), score(L, dense), tolerance = 1e-7)
    if (n < 7L) {
      expect_null(compact$Sigma)
      expect_null(compact$Sigma0)
    }
    expect_identical(score(L * 0, compact), -Inf)
    if (rank > 1L) {
      L[, 2L] <- L[, 1L]
      expect_equal(score(L, compact), score(L, dense), tolerance = 2e-6)
    }
  }
})

test_that("cached full-support loading roots retain the ridge on the nullspace", {
  set.seed(213)
  e <- private(".egcar_engine")(egcar_control(backend = "R"))
  for (n in c(6L, 35L)) {
    prep <- egcar_prepare(list(matrix(rnorm(n * 12L), n, 12L),
                              matrix(rnorm(n * 8L), n, 8L)))$prep
    for (ridge in c(1e-4, .01)) {
      selected <- seq_len(prep$p)
      factors <- e$egcar_loading_factors(prep, selected, ridge)
      scale_diag <- mean(unlist(lapply(prep$S_kk, diag)))
      for (k in seq_len(prep$K)) {
        ev <- eigen(e$symmetrize(prep$S_kk[[k]]), symmetric = TRUE)
        d <- pmax(ev$values + ridge * scale_diag, 1e-10)
        half <- tcrossprod(sweep(ev$vectors, 2L, sqrt(d), "*"), ev$vectors)
        inv <- tcrossprod(sweep(ev$vectors, 2L, 1 / sqrt(d), "*"), ev$vectors)
        expect_equal(e$egcar_apply_loading_factor(factors$blocks[[k]], diag(prep$p_list[[k]])), half, tolerance = 1e-8)
        expect_equal(e$egcar_apply_loading_factor(factors$blocks[[k]], diag(prep$p_list[[k]]), TRUE), inv, tolerance = 1e-8)
      }
    }
  }
})

test_that("compact and public group warm starts give the same path", {
  set.seed(214)
  p <- egcar_prepare(list(matrix(rnorm(12 * 17), 12, 17),
                          matrix(rnorm(12 * 8), 12, 8), matrix(rnorm(12 * 4), 12, 4)))$prep
  for (backend in c("cpp", "R", "reference")) {
    e <- private(".egcar_engine")(egcar_control(backend = backend))
    a <- e$fit_l21_admm(p, .03, max_iter = 30L, abs_tol = 0, rel_tol = 0,
                         compact_state = TRUE)
    b <- e$fit_l21_admm(p, .03, max_iter = 30L, abs_tol = 0, rel_tol = 0)
    expect_equal(a$C_hat, b$C_hat, tolerance = 1e-10)
    if (backend != "reference") {
      expect_null(a$G); expect_null(a$V)
      expect_true(all(c("Hk", "Hl", "a") %in% names(a)))
    }
    aa <- e$fit_l21_admm(p, .01, max_iter = 35L, abs_tol = 0, rel_tol = 0,
                         init = e$egcar_warm_state(a, TRUE))
    bb <- e$fit_l21_admm(p, .01, max_iter = 35L, abs_tol = 0, rel_tol = 0,
                         init = e$egcar_warm_state(b, TRUE))
    expect_equal(aa$C_hat, bb$C_hat, tolerance = 1e-9)
    expect_equal(aa$G, bb$G, tolerance = 1e-9)
    expect_equal(aa$V, bb$V, tolerance = 1e-9)
  }
})

test_that("CV task environments contain only fitting inputs", {
  ctl <- egcar_control(backend = "R")
  e <- private(".egcar_engine")(ctl)
  fun <- private(".egcar_cv_fold_task")(c(.03, .01), 1L, "l21", ctl, e)
  expect_setequal(ls(environment(fun)), c("lambda", "rank", "penalty", "control", "e", "order_path"))
  expect_identical(parent.env(environment(fun)), asNamespace("egcar"))
})

test_that("wide-fold CV preserves selection against dense-reference fitting", {
  set.seed(215)
  X <- list(matrix(rnorm(14 * 20), 14, 20), matrix(rnorm(14 * 10), 14, 10))
  folds <- egcar_cv_data(X, nfolds = 2L)
  for (penalty in c("l11", "l21")) {
    ctl <- egcar_control(max_iter = 50L, max_iter_cv = 40L,
      abs_tol = 0, rel_tol = 0, partial_eigen = FALSE)
    a <- egcar_cv(folds, 1L, penalty, c(.03, .01), control = ctl)
    ctl$backend <- "reference"
    b <- egcar_cv(folds, 1L, penalty, c(.03, .01), control = ctl)
    expect_equal(a$cv_fold_table$score, b$cv_fold_table$score, tolerance = 1e-6)
    expect_equal(a$lambda, b$lambda)
    expect_identical(a$fold_id, b$fold_id)
  }
})

test_that("loading-cache byte limits do not change fitted factors", {
  set.seed(216)
  prep <- egcar_prepare(list(matrix(rnorm(40 * 8), 40, 8), matrix(rnorm(40 * 7), 40, 7)))$prep
  ctl <- egcar_control(backend = "R", loading_cache_max_bytes = 0)
  uncached <- private(".egcar_engine")(ctl)
  a <- uncached$egcar_loading_factors(prep, seq_len(prep$p), 1e-4)
  expect_equal(length(prep$loading_factor_cache), 0L)
  ctl$loading_cache_max_bytes <- 100000
  cached <- private(".egcar_engine")(ctl)
  b <- cached$egcar_loading_factors(prep, seq_len(prep$p), 1e-4)
  expect_equal(a, b, tolerance = 0)
  expect_equal(length(prep$loading_factor_cache), 1L)
  entries <- as.list(prep$loading_factor_cache)
  expect_lte(sum(vapply(entries, function(x) as.numeric(object.size(x)), numeric(1L))), 100000)
})

test_that("compact fold payloads also work in multisession workers", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  set.seed(217)
  X <- list(matrix(rnorm(12 * 7), 12, 7), matrix(rnorm(12 * 5), 12, 5))
  folds <- egcar_cv_data(X, nfolds = 2L)
  ctl <- egcar_control(max_iter = 12L, max_iter_cv = 10L,
    abs_tol = 0, rel_tol = 0, partial_eigen = FALSE)
  old <- Sys.getenv("EGCAR_MULTICORE", unset = NA_character_)
  Sys.setenv(EGCAR_MULTICORE = "0")
  tryCatch({
    for (penalty in c("l11", "l21")) {
      a <- egcar_cv(folds, 1L, penalty, c(.03, .01), workers = 1L, control = ctl)
      b <- egcar_cv(folds, 1L, penalty, c(.03, .01), workers = 2L, control = ctl)
      expect_equal(a$cv_fold_table$score, b$cv_fold_table$score, tolerance = 1e-8)
      expect_equal(a$lambda, b$lambda)
      expect_identical(a$fold_id, b$fold_id)
      if (requireNamespace("parallelly", quietly = TRUE) && parallelly::supportsMulticore()) {
        Sys.setenv(EGCAR_MULTICORE = "1")
        forked <- egcar_cv(folds, 1L, penalty, c(.03, .01), workers = 2L, control = ctl)
        expect_equal(a$cv_fold_table$score, forked$cv_fold_table$score, tolerance = 1e-8)
        expect_equal(a$lambda, forked$lambda)
        expect_identical(a$fold_id, forked$fold_id)
        Sys.setenv(EGCAR_MULTICORE = "0")
      }
    }
  }, finally = {
    if (is.na(old)) Sys.unsetenv("EGCAR_MULTICORE") else Sys.setenv(EGCAR_MULTICORE = old)
  })
})

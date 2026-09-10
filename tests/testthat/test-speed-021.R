test_that("compiled EGCAR preserves full-space reference updates and inputs", {
  cases <- list(list(n = 40L, p = c(4L, 4L, 4L)),
                list(n = 35L, p = c(12L, 15L, 18L)),
                list(n = 12L, p = c(18L, 4L, 9L)))
  for (i in seq_along(cases)) {
    z <- cases[[i]]
    sim <- egcar_simulate(n = z$n, p_list = z$p, rank = 1L,
                          active_per_view = 2L, signal = 0.8, seed = 170L + i)
    if (i == 3L) sim$views[[2L]][, 1L] <- 0
    x <- egcar_prepare(sim$views)
    before <- serialize(x$prep$egcar_context, NULL)
    for (penalty in c("l11", "l21")) {
      ctl <- egcar_control(backend = "cpp", max_iter = 100L,
        abs_tol = 0, rel_tol = 0, check_every = 5L,
        partial_eigen = FALSE, keep_history = penalty == "l21")
      a <- egcar_fit(x, rank = 1L, penalty = penalty, lambda = 0.01, control = ctl)
      ctl$backend <- "R"
      b <- egcar_fit(x, rank = 1L, penalty = penalty, lambda = 0.01, control = ctl)
      fields <- if (penalty == "l11") c("C", "Z", "H") else c("C", "G", "V")
      for (nm in fields) {
        aa <- unlist(a$solver[[nm]], use.names = FALSE)
        bb <- unlist(b$solver[[nm]], use.names = FALSE)
        expect_equal(aa, bb, tolerance = 1e-8, info = paste(i, penalty, nm))
      }
      expect_identical(a$iterations, b$iterations)
      expect_identical(a$converged, b$converged)
      expect_identical(serialize(x$prep$egcar_context, NULL), before)
      init_before <- serialize(a$solver, NULL)
      ctl$backend <- "cpp"
      invisible(egcar_fit(x, 1L, penalty, 0.005, ctl, init = a))
      expect_identical(serialize(a$solver, NULL), init_before)
      expect_identical(serialize(x$prep$egcar_context, NULL), before)
    }
  }
})

test_that("compiled EGCAR CV keeps the same folds, candidates and held-out score", {
  sim <- egcar_simulate(n = 45L, p_list = c(5L, 5L, 5L),
                        rank = 1L, active_per_view = 2L, seed = 131L)
  shared <- egcar_cv_data(sim$views, nfolds = 3L, seed = 151L)
  for (penalty in c("l11", "l21")) {
    ctl <- egcar_control(backend = "cpp", max_iter = 150L, max_iter_cv = 100L,
                         abs_tol = 0, rel_tol = 0, partial_eigen = FALSE)
    a <- egcar_cv(shared, 1L, penalty, c(0.005, 0.02), workers = 1L, control = ctl)
    ctl$backend <- "R"
    b <- egcar_cv(shared, 1L, penalty, c(0.005, 0.02), workers = 1L, control = ctl)
    expect_identical(a$fold_id, shared$fold_id)
    expect_identical(a$fold_id, b$fold_id)
    expect_identical(a$cv_table$lambda, b$cv_table$lambda)
    expect_equal(a$cv_fold_table$score, b$cv_fold_table$score, tolerance = 1e-7)
    expect_identical(a$lambda, b$lambda)
    expect_error(egcar_cv(shared, penalty = penalty, lambda = c(0, 0.01)), "zero")
  }
})

test_that("the single-configuration all-method checker is installed and parseable", {
  path <- system.file("examples", "06_small_all_methods_check.R", package = "egcar")
  expect_true(nzchar(path))
  expect_error(parse(path), NA)
})

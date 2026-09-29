test_that("SGCA defaults use one six-hour budget and paper iteration defaults", {
  b <- benchmark_control()
  expect_equal(b$sgca_time_limit, 21600)
  expect_identical(b$sgca_init_max_iter, 1000L)
  expect_identical(b$sgca_tgd_max_iter, 15000L)
  c <- egcar_experiment_config()
  expect_equal(c$sgca_time_limit, 21600)
  expect_identical(c$sgca_max_iter_init, 1000L)
  expect_identical(c$sgca_max_iter_tgd, 15000L)
  for (value in list(0, -1, NA_real_, NaN, c(1, 2)))
    expect_error(benchmark_control(sgca_time_limit = value), "sgca_time_limit")
  expect_error(benchmark_control(sgca_tgd_max_iter = 1.5), "sgca_tgd_max_iter")
  expect_error(egcar_experiment_config(sgca_time_limit = 0), "sgca_time_limit")
  expect_equal(benchmark_control(sgca_time_limit = Inf)$sgca_time_limit, Inf)
  old <- unclass(c); old$sgca_time_limit <- NULL
  expect_equal(egcar:::validate_egcar_experiment_config(old)$sgca_time_limit, 21600)
})

test_that("unlimited TGD stops on its numerical criterion", {
  out <- egcar:::sgca_tgd_penalized(diag(c(3, 2, 1)), diag(3),
    matrix(c(1, 0, 0), 3, 1), rank = 1, k = 3, lambda = 1, max_iter = Inf, stopping = "absolute_change")
  expect_true(out$converged)
  expect_equal(out$iterations, 1)
  expect_equal(out$absolute_change, 0)
  expect_equal(abs(out$L[, 1]), c(1, 0, 0))
})

test_that("the supervisor distinguishes success, errors, and elapsed timeout", {
  run <- egcar:::.egcar_run_with_deadline
  good <- run(function() 42L, list(), 15)
  expect_false(good$timed_out)
  expect_equal(good$result, 42L)
  expect_error(run(function() stop("test solver failure"), list(), 15), "test solver failure")
  late <- run(function() { Sys.sleep(5); 42L }, list(), 0.2)
  expect_true(late$timed_out)
  expect_null(late$result)
  expect_gte(late$elapsed, 0.2)
  expect_lt(late$elapsed, 4)
})

test_that("supervised and in-process SGCA have the same numerical results", {
  set.seed(213)
  sim <- egcar_simulate(n = 24, p_list = c(3, 3, 3), active_per_view = 2)
  shared <- egcar_cv_data(sim$views, nfolds = 2, seed = 7)
  b <- benchmark_control(sgca_time_limit = Inf,
    sgca_init_max_iter = 5, sgca_tgd_max_iter = 10)
  run <- function(b, workers = 1L) sgca_cv(shared, rank = 1, k_grid = 6,
    rho_grid = .01, lambda_grid = 1, benchmarks = b, workers = workers)
  old <- run(b)
  b$sgca_time_limit <- 30
  new <- run(b)
  expect_false(new$timed_out)
  expect_identical(new$status, old$status)
  expect_identical(new$converged, old$converged)
  expect_equal(new$L, old$L, tolerance = 1e-12)
  expect_equal(new$cv_table, old$cv_table, tolerance = 1e-12)
  expect_gte(new$total_time, old$fit_time + old$tuning_time - 0.2)
  expect_equal(new$total_time, new$wall_time, tolerance = 1e-8)
  b$fast_sgca_initializer <- FALSE
  reference <- run(b)
  expect_equal(reference$L, new$L, tolerance = 1e-8)
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)) {
    b$fast_sgca_initializer <- TRUE
    parallel <- run(b, 2L)
    expect_equal(parallel$L, new$L, tolerance = 1e-10)
    expect_equal(parallel$cv_table, new$cv_table, tolerance = 1e-10)
  }
  b$sgca_time_limit <- 0.01
  timed <- expect_warning(run(b), NA)
  expect_equal(timed$status, "time_limit")
  expect_true(timed$timed_out)
  expect_false(timed$converged)
  expect_null(timed$L)
  expect_null(timed$best)
  expect_equal(timed$mean_loss, Inf)
  expect_match(timed$error, "did not converge in real time")
  expect_equal(timed$total_time, timed$wall_time, tolerance = 1e-8)
})

test_that("a timeout kills spawned R descendants as well as the coordinator", {
  skip_if_not(.Platform$OS.type == "unix" && dir.exists("/proc"))
  path <- tempfile("sgca-descendant-")
  on.exit(unlink(path), add = TRUE)
  late <- egcar:::.egcar_run_with_deadline(function(path) {
    p <- callr::r_bg(function() Sys.sleep(30), supervise = TRUE)
    writeLines(as.character(p$get_pid()), path)
    Sys.sleep(30)
  }, list(path = path), 2)
  expect_true(late$timed_out)
  expect_true(file.exists(path))
  pid <- as.integer(readLines(path))
  proc <- file.path("/proc", as.character(pid), "stat")
  alive <- file.exists(proc)
  if (alive) {
    state <- sub("^.*\\) ([A-Z]).*$", "\\1", readLines(proc, warn = FALSE))
    alive <- !state %in% c("Z", "X")
  }
  expect_false(alive)
})

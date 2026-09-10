test_that("only separate EGCAR families are exposed and CV excludes zero", {
  exports <- getNamespaceExports("egcar")
  expect_true(all(c("EGCAR_L11", "EGCAR_L21", "EGCAR_L11_CV", "EGCAR_L21_CV",
    "EGCAR_L11_Rate", "EGCAR_L21_Rate", "run_egcar_experiments") %in% exports))
  cfg <- egcar_experiment_config()
  expect_equal(cfg$rho_e_cv_grid, 10^seq(-5, 4))
  expect_equal(cfg$lambda_g_cv_grid, 10^seq(-5, 4))
  expect_error(egcar_experiment_config(rho_e_cv_grid = c(0, 0.01)), "strictly positive")
  expect_error(egcar_experiment_config(lambda_g_cv_grid = c(0, 0.01)), "strictly positive")
  sim <- egcar_simulate(n = 16, p_list = c(3, 3, 3), active_per_view = 2)
  expect_error(EGCAR_L11_CV(sim$views, lambda = c(0, 0.01)), "strictly positive")
  expect_error(EGCAR_L21_CV(sim$views, lambda = c(0, 0.01)), "strictly positive")
})

test_that("the SGCA initializer is bundled with all its required helpers", {
  ns <- asNamespace("egcar")
  expect_true(all(vapply(c("sgca_init_fixed", "updatePi", "updateH", "Soft"),
    exists, logical(1), envir = ns, inherits = FALSE)))
  e <- get(".egcar_engine", ns)(egcar_control(backend = "R"))
  A <- diag(5) + 0.05
  B <- diag(seq(1, 2, length.out = 5))
  reference <- e$get_sgca_initializer()
  z <- reference(A, B, rho = 0.05, K = 2, maxiter = 8, epsilon = 1e-12)
  prep <- e$sgca_prepare_initializer(A, B, reference)
  fast <- e$sgca_init_cached(prep, rho = 0.05, K = 2, maxiter = 8, epsilon = 1e-12)
  expect_equal(fast$Pi, z$Pi, tolerance = 1e-8)
  expect_equal(fast$H, z$H, tolerance = 1e-8)
  expect_equal(fast$Gamma, z$Gamma, tolerance = 1e-8)
  expect_equal(fast$iteration, z$iteration)
})

test_that("a tiny full experiment uses the configured methods and outputs", {
  cfg <- egcar_experiment_config(p_list = c(3L, 3L, 3L), n_grid = 18L, rank_grid = 1L,
    active_per_view = 2L, n_folds = 2L,
    rho_e_cv_grid = c(0.001, 0.01), lambda_g_cv_grid = c(0.001, 0.01),
    max_iter_cv = 10L, max_iter_final = 20L, oracle1_max_iter = 20L,
    run_external_benchmarks = FALSE, make_plots = FALSE, make_loading_plots = FALSE)
  outdir <- tempfile("egcar-test-")
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  ans <- suppressWarnings(run_egcar_experiments(outdir, config = cfg, backend = "R"))
  expected <- c("Oracle1-population", "Oracle2-support", "EGCAR-L11-rate", "EGCAR-L11-CV",
    "EGCAR-L21-rate", "EGCAR-L21-CV", "SGCA", "RGCCA", "SGCCA", "MultiCCA")
  expect_setequal(ans$results$method, expected)
  expect_true(file.exists(ans$checkpoint))
  expect_true(all(ans$cv_results$rho_e[ans$cv_results$method == "EGCAR-L11-CV"] > 0))
  expect_true(all(ans$cv_results$lambda_g[ans$cv_results$method == "EGCAR-L21-CV"] > 0))
  expect_equal(nrow(ans$egcar_cv_fold_results), 8L)
})

test_that("L21 defaults match the frozen pilot while explicit multipliers are preserved", {
  recorded <- as.numeric(readLines(system.file("validation", "l21_rate_0214",
                                               "selected_multiplier.txt", package = "egcar")))
  expect_length(recorded, 1L)
  cfg <- egcar_experiment_config()
  expect_equal(cfg$rate_c_g, recorded)
  expect_equal(cfg$rate_c_e, 1)
  expect_equal(egcar_experiment_config(rate_c_g = 1)$rate_c_g, 1)
  sim <- egcar_simulate(24, p_list = c(3, 4, 5), active_per_view = 2, seed = 214)
  prepared <- egcar_prepare(sim$views)
  a <- EGCAR_L21_Rate(prepared)
  b <- egcar_rate(prepared, penalty = "l21")
  explicit <- egcar_rate(prepared, penalty = "l21", multiplier = 1)
  base <- sqrt((9 + log(12)) / 24)
  direct <- egcar_fit(prepared, penalty = "l21", lambda = base)
  expect_equal(a$rate_multiplier, recorded)
  expect_identical(a$rate_multiplier_source, "independent_pilot_0214")
  expect_equal(a$lambda, recorded * base)
  expect_equal(a$C, b$C, tolerance = 1e-12)
  expect_equal(explicit$C, direct$C, tolerance = 1e-12)
  expect_identical(explicit$rate_multiplier_source, "explicit")
  l11 <- egcar_rate(prepared, penalty = "l11")
  expect_equal(l11$lambda, sqrt(log(12) / 24))
  expect_equal(l11$C, EGCAR_L11_Rate(prepared)$C, tolerance = 1e-12)
  expect_error(egcar_rate(prepared, penalty = "l21", multiplier = -1), "multiplier")
})

test_that("the zero certificate matches the analytic two-scalar optimum", {
  prep <- egcar_prepare(list(matrix(c(-1, 1), 2, 1), matrix(c(-1, 1), 2, 1)))
  expect_equal(egcar:::.egcar_l21_zero_bound(prep$prep), 0.5)
  control <- egcar_control(abs_tol = 1e-10, rel_tol = 1e-10)
  nonzero <- egcar_fit(prep, penalty = "l21", lambda = .49, control = control)
  zero <- egcar_fit(prep, penalty = "l21", lambda = .51, control = control)
  expect_equal(as.numeric(nonzero$C[[1L]]), .02, tolerance = 1e-7)
  expect_identical(nonzero$status, "ok")
  expect_equal(as.numeric(zero$C[[1L]]), 0, tolerance = 1e-9)
  expect_identical(zero$status, "invalid_loading")
  expect_true(zero$converged)
  rate <- egcar_rate(prep, penalty = "l21", multiplier = 10, control = control)
  expect_true(rate$rate_zero_certified)
  expect_equal(rate$rate_zero_bound, .5)
})

test_that("the tiled bound includes both endpoints and partial final tiles", {
  set.seed(214)
  X <- list(matrix(rnorm(6 * 257), 6, 257), matrix(rnorm(6 * 260), 6, 260),
            matrix(rnorm(6 * 3), 6, 3))
  prep <- egcar_prepare(X)$prep
  S <- matrix(0, 520, 520)
  ids <- egcar:::make_block_indices(prep$p_list)
  for (j in seq_len(nrow(prep$edge_table))) {
    e <- prep$edge_table[j, ]
    S[ids[[e$k]], ids[[e$l]]] <- prep$S_kl[[e$key]]
    S[ids[[e$l]], ids[[e$k]]] <- t(prep$S_kl[[e$key]])
  }
  expect_equal(egcar:::.egcar_l21_zero_bound(prep), .5 * max(sqrt(rowSums(S^2))), tolerance = 1e-12)
})

test_that("experiment status distinguishes solver convergence and loading validity", {
  sim <- egcar_simulate(18, p_list = c(3, 3, 3), rank = 2, active_per_view = 2, seed = 214)
  zero_C <- lapply(sim$population$Cstar, function(x) x * 0)
  evaluate <- function(loading, C = zero_C, converged = TRUE, status = "ok", error = NA_character_) {
    egcar:::evaluate_method("EGCAR-L21-rate", C, loading, sim$population,
      rank = 2, n = 18, rep_id = 1, lambda_g = 100,
      converged = converged, status = status, error_message = error, rate_zero_bound = .5)
  }
  z <- evaluate(list(valid = FALSE, reason = "fewer selected rows than rank"))
  expect_identical(z$status, "invalid_loading")
  expect_true(z$converged)
  expect_false(z$loading_valid)
  expect_true(z$zero_solution)
  expect_true(z$rate_zero_certified)
  expect_equal(z$selected_rows, 0)
  expect_true(is.na(z$subspace_sigma0))
  expect_match(z$error_message, "fewer selected rows")
  deficient <- evaluate(list(valid = TRUE, L = matrix(1, 9, 2)), C = NULL)
  expect_identical(deficient$status, "invalid_loading")
  expect_false(deficient$loading_valid)
  valid <- evaluate(list(valid = TRUE, L = sim$truth), C = NULL, converged = FALSE)
  expect_identical(valid$status, "not_converged")
  expect_true(valid$loading_valid)
  for (status in c("error", "time_limit", "skipped")) {
    failure <- evaluate(NULL, C = NULL, converged = FALSE, status = status, error = "original reason")
    expect_identical(failure$status, status)
    expect_identical(failure$error_message, "original reason")
  }
})

test_that("the experiment writes invalid L21-rate diagnostics to its failure CSV", {
  cfg <- egcar_experiment_config(p_list = c(3L, 3L, 3L), n_grid = 18L, rank_grid = 1L,
    active_per_view = 2L, n_folds = 2L, rho_e_cv_grid = .01, lambda_g_cv_grid = .01,
    rate_c_g = 1e5, max_iter_cv = 10L, max_iter_final = 50L, oracle1_max_iter = 20L,
    run_external_benchmarks = FALSE, make_plots = FALSE, make_loading_plots = FALSE,
    save_compact_loadings = FALSE)
  outdir <- tempfile("egcar-rate-")
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  ans <- suppressWarnings(run_egcar_experiments(outdir, config = cfg))
  row <- ans$results[ans$results$method == "EGCAR-L21-rate", ]
  expect_identical(row$status, "invalid_loading")
  expect_false(row$loading_valid)
  expect_true(row$zero_solution)
  expect_true(row$rate_zero_certified)
  failures <- read.csv(file.path(outdir, "benchmark_failures.csv"))
  expect_true("EGCAR-L21-rate" %in% failures$method)
})

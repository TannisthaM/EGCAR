test_that("legacy experiment configs receive retention defaults", {
  cfg <- egcar_experiment_config(
    p_list = c(4L, 4L, 4L), n_grid = 60L, rank_grid = 1L,
    active_per_view = 2L
  )
  legacy <- unclass(cfg)
  legacy$save_cv_fold_results <- NULL
  legacy$save_loading_data <- NULL
  legacy$save_compact_loadings <- NULL
  legacy$retain_benchmark_fits <- NULL
  validator <- getFromNamespace("validate_egcar_experiment_config", "egcar")
  got <- validator(legacy)
  expect_false(got$save_cv_fold_results)
  expect_false(got$save_loading_data)
  expect_true(got$save_compact_loadings)
  expect_false(got$retain_benchmark_fits)
})

test_that("memory-first experiment defaults are compact", {
  cfg <- egcar_experiment_config()
  expect_false(cfg$save_fits)
  expect_false(cfg$save_cv_fold_results)
  expect_false(cfg$save_loading_data)
  expect_true(cfg$save_compact_loadings)
  expect_false(cfg$retain_benchmark_fits)
  expect_false(cfg$make_loading_plots)
  expect_lte(cfg$loading_factor_cache_max, 4L)
})

test_that("memory retention flags can be restored explicitly", {
  cfg <- egcar_experiment_config(
    save_fits = TRUE,
    save_cv_fold_results = TRUE,
    save_loading_data = TRUE,
    save_compact_loadings = FALSE,
    retain_benchmark_fits = TRUE,
    loading_factor_cache_max = 128L
  )
  expect_true(cfg$save_fits)
  expect_true(cfg$save_cv_fold_results)
  expect_true(cfg$save_loading_data)
  expect_false(cfg$save_compact_loadings)
  expect_true(cfg$retain_benchmark_fits)
  expect_equal(cfg$loading_factor_cache_max, 128L)
})

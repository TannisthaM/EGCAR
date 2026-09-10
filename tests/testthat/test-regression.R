test_that("the updated numerical and public-API regression suite passes", {
  script <- system.file("validation", "core.R", package = "egcar")
  expect_true(nzchar(script))
  expect_error(sys.source(script, envir = new.env(parent = globalenv())), NA)
})

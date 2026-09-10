test_that("optional comparison and parallel regression checks pass", {
  skip_if_not(identical(Sys.getenv("EGCAR_TEST_OPTIONAL", "false"), "true"),
              "Set EGCAR_TEST_OPTIONAL=true to run optional comparisons.")
  script <- system.file("validation", "optional.R", package = "egcar")
  expect_true(nzchar(script))
  expect_error(sys.source(script, envir = new.env(parent = globalenv())), NA)
})

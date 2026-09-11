test_that("all original runnable examples remain installed", {
  examples <- c("00_install_dependencies.R", "01_quickstart.R", "02_all_methods_cv.R",
                "03_parallel_cv.R", "04_backend_validation.R", "05_backend_timing.R")
  for (name in examples) {
    path <- system.file("examples", name, package = "egcar")
    expect_true(nzchar(path), info = name)
    expect_error(parse(path), NA, info = name)
  }
})

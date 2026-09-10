test_that("every relocated numerical-engine function resolves in the namespace", {
  ns <- asNamespace("egcar")
  registry <- get(".egcar_engine_functions", envir = ns, inherits = FALSE)
  expect_false(anyDuplicated(registry) > 0L)
  expect_true(all(vapply(registry, function(nm) {
    exists(nm, envir = ns, inherits = FALSE) &&
      is.function(get(nm, envir = ns, inherits = FALSE))
  }, logical(1))))
  e <- get(".egcar_engine", envir = ns)(egcar_control(backend = "R"))
  expect_true(all(vapply(registry, function(nm) {
    identical(environment(get(nm, envir = e, inherits = FALSE)), e)
  }, logical(1))))
})

test_that("all fitting and comparison entry points remain exported", {
  expected <- c("egcar_fit", "egcar_rate", "egcar_path", "egcar_cv",
                "egcar_cv_data", "sgca_cv", "rgcca_cv", "sgcca_cv", "multicca_cv")
  expect_true(all(expected %in% getNamespaceExports("egcar")))
  expect_true(all(vapply(expected, function(nm) is.function(getExportedValue("egcar", nm)),
                         logical(1))))
})

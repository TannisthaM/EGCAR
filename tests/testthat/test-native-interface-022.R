test_that("solver state dimensions are checked, not silently reconstructed", {
  check <- getFromNamespace(".egcar_check_solver_state", "egcar")
  context <- list(p_list = c(2L, 3L), edge_k = 1L, edge_l = 2L)
  state <- list(C = list(matrix(1:6, 2L, 3L)),
                Z = list(matrix(0, 2L, 3L)), H = list(matrix(0, 2L, 3L)))
  before <- serialize(state, NULL)
  expect_true(check(state, context, FALSE, "cpp"))
  expect_identical(serialize(state, NULL), before)
  bad <- state; bad$C[[1L]] <- as.numeric(bad$C[[1L]])
  expect_error(check(bad, context, FALSE, "cpp"), "state\\$C\\[\\[1\\]\\].*2 x 3")
  bad <- state; bad$C[[1L]] <- t(bad$C[[1L]])
  expect_error(check(bad, context, FALSE, "cpp"), "3 x 2")
  bad <- state; bad$C[[1L]] <- array(1:6, c(2L, 3L, 1L))
  expect_error(check(bad, context, FALSE, "cpp"), "2 x 3 x 1")
  bad <- state; bad$Z <- numeric(6L)
  expect_error(check(bad, context, FALSE, "cpp"), "state\\$Z must contain")
  bad <- state; bad$H <- NULL
  expect_error(check(bad, context, FALSE, "cpp"), "required named state lists")
})

test_that("a mismatched DLL fails before row norms or loading extraction", {
  e <- getFromNamespace(".egcar_engine", "egcar")(egcar_control(backend = "cpp"))
  sim <- egcar_simulate(n = 12L, p_list = c(2L, 3L), active_per_view = 1L, toeplitz_rho = 0.5)
  prep <- egcar_prepare(sim$views)$prep
  # Override only the test engine's native function; never alter the namespace.
  e$egcar_native_solve <- function(context, state, controls, group, verbose) list(state = state)
  ctl <- e$egcar_controls(0.01, 1, 2L, 0, 0, FALSE, 10, 2, 10L, 1L, FALSE)
  expect_error(e$egcar_run_solver(prep, ctl, NULL, FALSE), "native matrix API mismatch")
})

test_that("the native matrix interface and signal-0.8 rank-1/2/5 paths work", {
  path <- system.file("validation", "matrix_interface.R", package = "egcar", mustWork = TRUE)
  checks <- new.env(parent = globalenv())
  sys.source(path, envir = checks)
  report <- checks$run_egcar_matrix_interface_check(workers = 1L, check_cv = TRUE, verbose = FALSE)
  expect_equal(nrow(report$native), 16L)
  expect_equal(nrow(report$fits), 30L)
  expect_equal(nrow(report$cv), 6L)
  expect_true(all(report$native$passed) && all(report$fits$passed) && all(report$cv$passed))
  expect_setequal(report$fits$rank, c(1L, 2L, 5L))
  expect_true(all(report$fits$signal == 0.8))
})

test_that("the focused matrix-interface example is installed and parseable", {
  path <- system.file("examples", "07_matrix_interface_check.R", package = "egcar", mustWork = TRUE)
  expect_error(parse(path), NA)
})

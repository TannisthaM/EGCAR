test_that("legacy dimensionless edge blocks are restored safely", {
  z <- list(
    p_list = c(2L, 3L, 2L),
    edge_k = c(1L, 1L, 2L),
    edge_l = c(2L, 3L, 3L),
    cols_k = list(1:3, 4:5, 1:2),
    cols_l = list(1:2, 1:2, 3:5)
  )
  mats <- list(matrix(seq_len(6), 2, 3), matrix(seq_len(4), 2, 2), matrix(seq_len(6), 3, 2))
  flat <- lapply(mats, as.numeric)
  n1 <- egcar_group_norms(z, mats)
  n2 <- egcar_group_norms(z, flat)
  expect_equal(n2, n1)

  state <- list(C = flat, Z = flat, H = flat)
  repaired <- .egcar_normalize_solver_state(state, z, FALSE)
  expect_true(all(vapply(repaired$C, is.matrix, logical(1))))
  expect_identical(dim(repaired$C[[1L]]), c(2L, 3L))
  expect_identical(dim(repaired$C[[2L]]), c(2L, 2L))
  expect_identical(dim(repaired$C[[3L]]), c(3L, 2L))
})

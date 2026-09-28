# Frozen by tools/calibrate_l21_rate_0214.R from 36 independent pilot cases.
# See inst/validation/l21_rate_0214/pilot_summary.csv and the calibration report.
# This is the best tested held-out-score coefficient, not a sparsity guarantee.
.egcar_l21_rate_multiplier <- 0.05

.egcar_default_rate_multiplier <- function(penalty) {
  if (identical(penalty, "l21")) .egcar_l21_rate_multiplier else 1
}

# Sufficient zero-solution certificate for the overlapping endpoint groups.
# Each cross-covariance coefficient occurs in two incident rows. Cauchy-Schwarz
# gives F(C)-F(0) >= (lambda - b) * sum_k ||M_k(C)||_2,1, with b below.
# lambda > b certifies the unique zero optimum. lambda < b is inconclusive.
# Tile existing edges instead of assembling the full covariance or squaring
# an entire large edge at once. Additional storage is O(p + 256^2).
.egcar_l21_zero_bound <- function(prep) {
  norms2 <- lapply(prep$p_list, numeric)
  for (e in seq_len(nrow(prep$edge_table))) {
    k <- prep$edge_table$k[[e]]
    l <- prep$edge_table$l[[e]]
    S <- prep$S_kl[[prep$edge_table$key[[e]]]]
    for (first_row in seq.int(1L, nrow(S), by = 256L)) {
      rows <- seq.int(first_row, min(first_row + 255L, nrow(S)))
      for (first_col in seq.int(1L, ncol(S), by = 256L)) {
        cols <- seq.int(first_col, min(first_col + 255L, ncol(S)))
        tile <- S[rows, cols, drop = FALSE]^2
        norms2[[k]][rows] <- norms2[[k]][rows] + rowSums(tile)
        norms2[[l]][cols] <- norms2[[l]][cols] + colSums(tile)
      }
    }
  }
  0.5 * sqrt(max(vapply(norms2, max, numeric(1L))))
}

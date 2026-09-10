library(egcar)
set.seed(104)
views <- list(matrix(rnorm(30 * 40), 30, 40), matrix(rnorm(30 * 12), 30, 12),
              matrix(rnorm(30 * 8), 30, 8))
# Include a constant column to exercise a singular covariance.
views[[1]][, 1] <- 0
prepared <- egcar_prepare(views)
for (penalty in c("l11", "l21")) {
  fits <- lapply(c("cpp", "R", "reference"), function(backend) {
    egcar_fit(prepared, rank = 1, penalty = penalty, lambda = 0.01,
      control = egcar_control(backend = backend, max_iter = 200,
        abs_tol = 0, rel_tol = 0, adaptive_mu = TRUE, check_every = 5))
  })
  baseline <- unlist(fits[[3]]$C, use.names = FALSE)
  for (j in 1:2) {
    difference <- max(abs(unlist(fits[[j]]$C, use.names = FALSE) - baseline))
    cat(penalty, c("cpp", "R")[[j]], "max operator difference:", difference, "\n")
    stopifnot(difference < 1e-6 * max(1, max(abs(baseline))))
  }
}

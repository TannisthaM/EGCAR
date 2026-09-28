# Small installed-package check; these are diagnostics, not timing claims.
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.11")
source(system.file("validation", "core.R", package = "egcar", mustWork = TRUE))
set.seed(211)
X <- list(matrix(rnorm(24 * 35), 24, 35), matrix(rnorm(24 * 8), 24, 8))
shared <- egcar_cv_data(X, nfolds = 3L)
ctl <- egcar_control(max_iter = 60L, max_iter_cv = 40L,
                      abs_tol = 0, rel_tol = 0, partial_eigen = FALSE)
for (penalty in c("l11", "l21")) {
  ctl$backend <- "cpp"
  a <- egcar_cv(shared, 1L, penalty, c(.03, .01), control = ctl)
  ctl$backend <- "R"
  b <- egcar_cv(shared, 1L, penalty, c(.03, .01), control = ctl)
  stopifnot(isTRUE(all.equal(a$cv_fold_table$score, b$cv_fold_table$score, tolerance = 1e-6)),
            identical(a$lambda, b$lambda))
}
cat("0.2.11 installed-package performance regression checks passed.\n")

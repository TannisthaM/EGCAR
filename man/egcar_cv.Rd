\name{egcar_cv}
\alias{egcar_cv}
\alias{print.egcar_cv}
\title{Cross-validate either accelerated EGCAR penalty}
\description{Evaluate a full one-dimensional direct-coefficient grid on common folds, select the highest average held-out generalized Rayleigh score, and refit on all observations.}
\usage{
egcar_cv(x, rank = 1L, penalty = c("l11", "l21"),
                     lambda = 10^seq(-5, 4), fold_id = NULL,
                     nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control())
}
\arguments{
\item{x}{Raw views or a reusable \code{egcar_cv_data} object.}
\item{rank}{Requested positive loading rank.}
\item{penalty}{Separate entrywise \code{"l11"} or row-group \code{"l21"} estimator.}
\item{lambda}{Direct coefficient grid, strictly positive; zero is excluded. No fold-size scaling is applied.}
\item{fold_id}{Optional supplied fold labels, or NULL to use labels already stored in x.}
\item{nfolds}{Number of folds when creating them from raw views.}
\item{seed}{Nonnegative integer seed; the caller's RNG state is restored.}
\item{workers}{Fold-level workers, capped at the number of folds. Values above one require future and future.apply. The previous future plan is restored.}
\item{control}{An \code{egcar_control} object or named argument list.}
}
\details{Candidates are visited in descending coefficient order within each fold
and use same-family warm starts. By default every fold must converge and have
a finite requested-rank score. Set \code{cv_require_convergence=FALSE} in
\code{egcar_control} only to reproduce the earlier finite-score eligibility rule.
Every fold must have a finite score for a candidate to be eligible. Exact score
ties favor stronger regularization. Zero is excluded from the CV grid for both EGCAR penalties.

The common loss is the negative generalized Rayleigh score, with the supplied
relative ridge and eigenvalue floor. It is not reconstruction MSE. Validation
covariances use training-mean centering only. If every candidate fails, the
result has status \code{"no_valid_cv"}, no selected penalty and no refit. This
explicit failure replaces the earlier script's fallback to an arbitrary
coefficient of one. Diagnostics are in \code{cv_fold_table}.

Preparation is reported separately and can be amortized over methods. Tuning
time includes fold fits and scoring, and may include worker startup. EGCAR's
final \code{fit_time} includes optimization through loading extraction.
\code{solver_time} and \code{loading_time} are components of that interval;
\code{total_time = tuning_time + fit_time} does not add loading twice.
\code{wall_time} measures the actual current public call, including preparation when performed. The time of loading
or compiling the installed package is not included.}
\value{An \code{egcar_cv} object containing \code{fit}, \code{L}, \code{best}, \code{lambda}, candidate and fold tables, supplied fold labels, status, controls and timings.}
\examples{
sim <- egcar_simulate(n = 24, p_list = c(3, 4, 5), active_per_view = 2)
shared <- egcar_cv_data(sim$views, nfolds = 2)
cv <- egcar_cv(shared, penalty = "l11", lambda = c(0.01, 0.1))
print(cv)
}

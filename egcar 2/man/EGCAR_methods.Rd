\name{EGCAR_L11}
\alias{EGCAR_L11}
\alias{EGCAR_L21}
\alias{EGCAR_L11_CV}
\alias{EGCAR_L21_CV}
\alias{EGCAR_L11_Rate}
\alias{EGCAR_L21_Rate}
\title{Explicit Single-Penalty EGCAR Estimators}
\description{Convenience entry points for separate entrywise and concatenated-row-group EGCAR fitting, positive-grid cross-validation, and rate-scaled fitting.}
\usage{
EGCAR_L11(x, rank = 1L, lambda = 0.01, control = egcar_control(), init = NULL)
EGCAR_L21(x, rank = 1L, lambda = 0.01, control = egcar_control(), init = NULL)
EGCAR_L11_CV(x, rank = 1L, lambda = 10^seq(-5, 4), fold_id = NULL,
  nfolds = 5L, seed = 1L, workers = 1L, control = egcar_control())
EGCAR_L21_CV(x, rank = 1L, lambda = 10^seq(-5, 4), fold_id = NULL,
  nfolds = 5L, seed = 1L, workers = 1L, control = egcar_control())
EGCAR_L11_Rate(x, rank = 1L, multiplier = 1, control = egcar_control())
EGCAR_L21_Rate(x, rank = 1L, multiplier = 1, control = egcar_control())
}
\arguments{
\item{x}{Raw views, a prepared object, or shared CV data.}
\item{rank}{Requested number of leading positive loading directions.}
\item{lambda}{A nonnegative fixed-fit coefficient; for CV, a vector of strictly positive direct coefficients.}
\item{control}{Controls from \code{egcar_control()}.}
\item{init}{Optional same-family warm start for a fixed fit.}
\item{fold_id}{Optional shared fold labels.}
\item{nfolds}{Number of folds when no labels are supplied.}
\item{seed}{Random seed for folds and reproducible CV.}
\item{workers}{Fold-level worker count.}
\item{multiplier}{Nonnegative multiplier of the inherited rate scaling.}
}
\details{L11 adds only the entrywise penalty. L21 adds only
\eqn{\lambda\sum_k\|M_k(C)\|_{2,1}}, with complete concatenated rows across all
incident edges. No combined-penalty or tied-CV estimator is provided.
CV rejects zero for both penalties. An unused penalty coefficient is zero by
construction; that is not an unregularized candidate. Zero-coefficient fixed
fits and the population oracle remain available. See \code{egcar_fit},
\code{egcar_cv}, and \code{egcar_rate} for numerical controls and return fields.}
\value{An \code{egcar_fit} or \code{egcar_cv} object, as appropriate.}
\examples{
sim <- egcar_simulate(n = 24, p_list = c(3, 4, 5), active_per_view = 2)
a <- EGCAR_L21(sim$views, lambda = 0.01)
print(a)
}

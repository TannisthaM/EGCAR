\name{egcar_prepare}
\alias{egcar_prepare}
\alias{egcar_folds}
\alias{egcar_cv_data}
\alias{print.egcar_prepared}
\alias{print.egcar_cv_data}
\title{Prepare views and reusable common cross-validation folds}
\description{Compute cached covariance blocks, spectral reductions, training means and, when requested, common folds for all methods.}
\usage{
egcar_prepare(views)

egcar_folds(n, nfolds = 5L, seed = 1L)

egcar_cv_data(views, fold_id = NULL, nfolds = 5L, seed = 1L)
}
\arguments{
\item{views}{A list of at least two finite numeric matrices with matching sample rows. Numeric data frames are accepted. Variable and block scaling are not performed. Missing observations and missing views are not supported.}
\item{n}{Number of observations.}
\item{nfolds}{Number of nonempty folds.}
\item{seed}{Nonnegative integer seed. The caller's RNG state is restored.}
\item{fold_id}{One nonmissing fold label per row. Supplied labels are never regenerated.}
}
\details{Every covariance uses divisor n for that split. Validation observations are centered at training means, not their own means. All training splits must have at least two observations. A prepared object stores covariance information and means but not the full raw data. A CV-data object retains raw data and full and fold preparations, so it can be reused for both EGCAR penalties and all comparison methods. Stored means and raw data may be sensitive; save these objects accordingly.}
\value{\code{egcar_prepare} returns an \code{egcar_prepared} object. \code{egcar_folds} returns integer fold labels. \code{egcar_cv_data} returns an \code{egcar_cv_data} object.}
\examples{
sim <- egcar_simulate(n = 30, p_list = c(4, 5, 6), active_per_view = 3)
shared <- egcar_cv_data(sim$views, nfolds = 3)
print(shared)
}

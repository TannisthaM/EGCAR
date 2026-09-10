\name{egcar_fit}
\alias{egcar_fit}
\alias{egcar_rate}
\alias{egcar_path}
\alias{print.egcar_fit}
\title{Fit accelerated entrywise-only or group-only EGCAR}
\description{Fit the pairwise regression loss with one of two separate sparsity penalties and obtain localized generalized-correlation loadings.}
\usage{
egcar_fit(x, rank = 1L, penalty = c("l11", "l21"), lambda = 0.01,
                      control = egcar_control(), init = NULL)

egcar_rate(x, rank = 1L, penalty = c("l11", "l21"), multiplier = 1,
                       control = egcar_control())

egcar_path(x, rank = 1L, penalty = c("l11", "l21"),
                       lambda = 10^seq(-5, 4), control = egcar_control())
}
\arguments{
\item{x}{Raw views, an \code{egcar_prepared} object, or an \code{egcar_cv_data} object (using its full data).}
\item{rank}{Number of leading positive loading directions requested.}
\item{penalty}{\code{"l11"} for entrywise sparsity or \code{"l21"} for the view-wise concatenated row-group penalty. They are not combined.}
\item{lambda}{A direct nonnegative penalty coefficient, or the coefficient vector for a path.}
\item{control}{An \code{egcar_control} object or a named subset of its arguments.}
\item{init}{An optional same-family fitted object or dimension-compatible ADMM-state list. Used as a warm start.}
\item{multiplier}{Nonnegative multiplier applied to the selected rate scale.}
}
\details{For each unordered view pair the smooth loss is
\deqn{\tfrac12\mathrm{tr}(C_{kl}^{T}S_{kk}C_{kl}S_{ll})-\langle S_{kl},C_{kl}\rangle.}
The entrywise version adds \eqn{\lambda\sum_{k<l}\|C_{kl}\|_{1,1}}.
The group version adds \eqn{\lambda\sum_k\|M_k(C)\|_{2,1}}, where each row
of \eqn{M_k(C)} concatenates all coefficients incident to one variable in view k.
The latter uses only C, G and V; it contains no entrywise proximal update.
The reduced-basis linear solve retains the full null-space contribution.
Loadings use the largest algebraic eigenvalues of the localized normalized operator.
If fewer than the requested number of positive directions exist, the operator fit
is returned with status \code{"invalid_loading"}; no lower rank is silently substituted.
The rate scales are \eqn{\sqrt{\log(p)/n}} and
\eqn{\sqrt{(d_{max}+\log(p))/n}}, respectively, where
\eqn{d_{max}=\max_k(p-p_k)}. The group scale is inherited as a numerical
benchmark normalization, not asserted to be an optimal theoretical rate.
Paths visit coefficients from strongest to weakest, return them in input order,
and reuse primal/scaled-dual states with the configured initial augmentation
parameter, matching the supplied script's convention.}
\value{\code{egcar_fit} and \code{egcar_rate} return \code{egcar_fit} objects containing stacked loadings \code{L}, view loadings, edge blocks \code{C}, the raw solver state, residuals, convergence status, selected rows, training means, controls and separated timings. A path returns the coefficient vector and corresponding fitted objects.}
\examples{
sim <- egcar_simulate(n = 30, p_list = c(3, 4, 5), active_per_view = 2)
fit <- egcar_fit(sim$views, penalty = "l21", lambda = 0.01)
print(fit)
}

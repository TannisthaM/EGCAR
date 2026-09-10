\name{egcar_simulation}
\alias{egcar_simulate}
\alias{egcar_distance}
\alias{egcar_oracle_support}
\alias{egcar_oracle_population}
\title{Simulation model, subspace distance and oracle diagnostics}
\description{Reproduce the supplied positive-semidefinite multiview latent-factor model and its two truth-informed oracle benchmarks.}
\usage{
egcar_simulate(n = 120L, p_list = c(15L, 15L, 15L), rank = 1L,
                           active_per_view = 5L, toeplitz_rho = c(0.5, 0.7, 0.9),
                           signal = 0.8, seed = 1L)

egcar_distance(A, B, rank = ncol(A))

egcar_oracle_support(x, active_local, rank = 1L, ridge = 1e-8,
                                 control = egcar_control())

egcar_oracle_population(population, rank = population$rank,
                                    control = egcar_control())
}
\arguments{
\item{n}{Number of observations.}
\item{p_list}{Positive view dimensions.}
\item{rank}{Latent rank or requested loading rank.}
\item{active_per_view}{Active variable count per view, restricted by each view dimension.}
\item{toeplitz_rho}{One Toeplitz correlation or one per view, all with absolute value below one.}
\item{signal}{Shared signal strictly between zero and one.}
\item{seed}{Nonnegative integer seed; caller RNG state is restored.}
\item{A, B}{Finite matrices with equal row counts; their column spaces are compared.}
\item{x}{Raw views, prepared data or shared CV data for the support oracle.}
\item{active_local}{List of known true local column supports. Must not be used as information in ordinary CV.}
\item{ridge}{Small support-oracle linear-solve ridge.}
\item{control}{An \code{egcar_control} object or named subset.}
\item{population}{Population component returned by \code{egcar_simulate}.}
}
\details{The model, covariance divisor, normalization and oracle definitions are
inherited from the supplied script. Oracle1 runs the original zero-penalty consensus
ADMM on exact population covariance blocks with both penalty coefficients equal
to zero; this diagnostic is intentionally not substituted with the accelerated
group-only splitting. Oracle2 solves on empirical covariance blocks restricted
to known true supports. Neither oracle performs cross-validation.
\code{egcar_distance} is the Frobenius sine-theta distance after Euclidean
orthonormalization, not an aligned coefficient difference. It returns NA for
numerically deficient rank. Premultiply both bases by the population metric
square root for the supplied Sigma0-subspace diagnostic.}
\value{Simulation returns \code{views}, \code{population}, \code{truth} and support information. Distance returns a scalar. Oracles return fitted objects.}
\examples{
sim <- egcar_simulate(n = 30)
oracle <- egcar_oracle_support(sim$views, sim$active_local)
egcar_distance(oracle$L, sim$truth)
}

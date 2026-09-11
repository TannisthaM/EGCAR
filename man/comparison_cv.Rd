\name{comparison_cv}
\alias{sgca_cv}
\alias{rgcca_cv}
\alias{sgcca_cv}
\alias{multicca_cv}
\title{Common-loss CV for SGCA, RGCCA, SGCCA and MultiCCA}
\description{Optional comparison wrappers retaining the supplied script's fitting algorithms, common validation score and fold-level parallelism.}
\usage{
sgca_cv(x, rank = 1L, k_grid = NULL,
                    rho_grid = c(0, 1e-3, 1e-2, 0.1, 0.5, 1),
                    lambda_grid = 10^seq(-5, 4), fold_id = NULL, nfolds = 5L,
                    seed = 1L, workers = 1L, control = egcar_control(),
                    benchmarks = benchmark_control())

rgcca_cv(x, rank = 1L, tau_grid = c(1e-6, 1e-3, 0.1, 0.25, 0.5, 0.75, 1),
                     fold_id = NULL, nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control(), benchmarks = benchmark_control())

sgcca_cv(x, rank = 1L, sparsity_grid = NULL, fold_id = NULL,
                     nfolds = 5L, seed = 1L, workers = 1L,
                     control = egcar_control(), benchmarks = benchmark_control())

multicca_cv(x, rank = 1L, penalty_grid = NULL, fold_id = NULL,
                        nfolds = 5L, seed = 1L, workers = 1L,
                        control = egcar_control(), benchmarks = benchmark_control())
}
\arguments{
\item{x}{Raw views or the same \code{egcar_cv_data} object used for EGCAR.}
\item{rank}{Target loading rank.}
\item{k_grid}{SGCA total row-support budget across all views. NULL uses valid members of 5, 10, 15, 20, 30 and p.}
\item{rho_grid}{SGCA direct initializer coefficients; not rate multipliers.}
\item{lambda_grid}{SGCA strictly positive penalized-gradient coefficients.}
\item{tau_grid}{RGCCA tied tau grid in [0,1].}
\item{sparsity_grid}{SGCCA tied sparsity bounds between the maximum of 1/sqrt(p_k) and one. NULL uses a dimension-valid grid.}
\item{penalty_grid}{MultiCCA tied L1 bounds between one and the minimum of sqrt(p_k). NULL uses a dimension-valid grid.}
\item{fold_id, nfolds, seed, workers}{Common-fold and parallel controls as in \code{egcar_cv}.}
\item{control}{EGCAR infrastructure controls, including scoped BLAS threads. EGCAR ADMM controls do not change comparison solvers.}
\item{benchmarks}{A \code{benchmark_control} object or named argument list.}
}
\details{SGCA uses the bundled initializer
\code{sgca_init_fixed}, from the TannisthaM/SGCA repository. Its initializer is
cached by fold and direct rho. The subsequent penalized thresholded-gradient
implementation is the one in the supplied script, including backtracking and
no per-iteration metric normalization; the supplied local-experiment CV and TGD routines are retained.
It tunes k, rho and lambda jointly.

RGCCA and SGCCA use the public \code{RGCCA::rgcca} solver with variable/block
scaling disabled. They use \code{astar} loadings mapping original centered data
to deflated components. MultiCCA uses the supplied cached-Gram implementation
and the installed PMA threshold helpers, or \code{PMA::MultiCCA} when selected.
External block signs are resolved using training observations only.

All wrappers minimize exactly the same held-out loss as EGCAR. The wrappers do
not use the comparison packages' own cross-validation criteria. A candidate
must have finite loss on every fold. A returned finite approximate fit does
not imply convergence; diagnostic flags are retained. These methods return
loading estimates, not an EGCAR cross-view operator.

Optional packages are never automatically installed. Run the included
\code{00_install_dependencies.R --benchmarks} example explicitly when needed.
For external methods loading extraction is included in \code{fit_time};
\code{loading_time} is therefore zero, not an additional omitted cost.}
\value{An \code{egcar_cv} object with \code{L}, \code{best}, the external \code{fit_full}, candidate and fold diagnostics, status, means, controls, labels and timings.}
\examples{
\dontrun{
sim <- egcar_simulate(n = 60)
shared <- egcar_cv_data(sim$views, nfolds = 3)
rgcca_cv(shared, tau_grid = c(0.1, 1))
sgcca_cv(shared, sparsity_grid = c(0.5, 1))
multicca_cv(shared, penalty_grid = c(2, 3))
sgca_cv(shared, k_grid = c(10, 15), rho_grid = 0.01, lambda_grid = 0.1)
}
}

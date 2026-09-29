\name{comparison_cv}
\alias{sgca_cv}
\alias{rgcca_cv}
\alias{sgcca_cv}
\alias{multicca_cv}
\title{Cross-validation for SGCA, RGCCA, SGCCA and MultiCCA}
\description{Paper-based SGCA and optional comparison wrappers using shared folds and fold-level parallelism.}
\usage{
sgca_cv(x, rank = 1L, k_grid = NULL,
                    rho_grid = NULL,
                    lambda_grid = 0.01, fold_id = NULL, nfolds = 5L,
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
\item{k_grid}{SGCA total row-support budget across all views. NULL uses valid members of 5, 10, ..., 100; p is used if none are valid.}
\item{rho_grid}{NULL uses the paper rate per training sample; numeric values are direct initializer coefficients.}
\item{lambda_grid}{SGCA strictly positive penalized-gradient coefficients.}
\item{tau_grid}{RGCCA tied tau grid in [0,1].}
\item{sparsity_grid}{SGCCA tied sparsity bounds between the maximum of 1/sqrt(p_k) and one. NULL uses a dimension-valid grid.}
\item{penalty_grid}{MultiCCA tied L1 bounds between one and the minimum of sqrt(p_k). NULL uses a dimension-valid grid.}
\item{fold_id, nfolds, seed, workers}{Common-fold and parallel controls as in \code{egcar_cv}.}
\item{control}{EGCAR infrastructure controls, including scoped BLAS threads. EGCAR ADMM controls do not change comparison solvers.}
\item{benchmarks}{A \code{benchmark_control} object or named argument list.}
}
\details{SGCA defaults follow Gao and Ma (2023), Algorithm 1 and Sections 5.1--5.2:
fixed 15000 thresholded-gradient updates, eta=0.001, lambda=0.01, a generalized
Fantope initializer with trace exactly rank, and rho=0.5*sqrt(log(p)/n_train).
Covariances use divisor n. Initial loading extraction uses leading algebraic
eigenvectors, followed by hard row thresholding, metric normalization and
rescaling. Each gradient update uses 2*eta and has no backtracking or intermediate
normalization. The final loading is metric-normalized.

The initializer uses the authors' R-code tolerance 0.005 and an enforced
1000-outer-iteration cap; these are implementation controls, not numerical
accuracy guarantees specified by the paper. The original R counter bug is not
reproduced. Initialization is cached within each fold and rho value.

By default only sparsity is tuned, on valid members of 5,10,...,100, using the
paper's held-out covariance trace score. On dimensions below the first valid
grid point, the public NULL-grid default uses all p rows. NULL rho_grid applies
the rate separately on each training fold and full refit; explicit numeric
rho_grid values retain their old DIRECT-coefficient meaning. A supplied
lambda_grid also overrides the paper default. Five folds are the default;
precomputed folds are respected.

benchmark_control(sgca_cv_score="common_loss") selects EGCAR's common score
instead; this is an explicit comparison variant. sgca_stopping="absolute_change"
allows early stopping at an absolute Frobenius iterate change of 1e-6 by default.
That variant still uses the paper's normalization and 2*eta update, so it is not
a byte-for-byte reproduction of the authors' R implementation.

For fixed iterations, status="ok" and completed=TRUE mean a finite full-rank
benchmark output completed the requested steps. converged separately requires
both numerical change tests. A finite capped initializer is allowed by the
benchmark default and its cap is reported. Use sgca_cv_require_convergence=TRUE
to reject such candidates. EGCAR's cv_require_convergence does not control SGCA.
Fold diagnostics and final diagnostics record actual rho, stage iterations,
change magnitudes, convergence and stop reasons, even when dense fits are dropped.

The separate six-hour budget covers the entire SGCA CV grid and refit, including
initialization, not each fit. A timeout returns status="time_limit" and no loading.
Shared data preparation precedes the budget. Finite budgets use a supervised R
process; sgca_time_limit=Inf disables supervision. Numerical failures are not
repaired by changing rank, adding an unrequested ridge or switching algorithms.

RGCCA and SGCCA use the public \code{RGCCA::rgcca} solver with variable/block
scaling disabled. They use \code{astar} loadings mapping original centered data
to deflated components. MultiCCA uses the supplied cached-Gram implementation
and the installed PMA threshold helpers, or \code{PMA::MultiCCA} when selected.
External block signs are resolved using training observations only.

RGCCA, SGCCA and MultiCCA use the common EGCAR validation loss. SGCA uses its selected score. A candidate
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
sgca_cv(shared)
}
}

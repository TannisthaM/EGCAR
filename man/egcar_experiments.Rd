\name{run_egcar_experiments}
\alias{run_egcar_experiments}
\alias{egcar_experiment_config}
\title{Run the Complete Local EGCAR Comparison Study}
\description{Reproduces the supplied local simulation design using installed package functions, positive EGCAR CV grids, and the original comparison-method CV routines.}
\usage{
egcar_experiment_config(p_list = c(15L, 15L, 15L),
  n_grid = c(30L, 45L, 100L, 1000L, 5000L, 10000L),
  rank_grid = c(1L, 2L, 5L), ...)
run_egcar_experiments(output_dir = "egcar_package_outputs",
  workers = 1L, n_reps = 1L, config = egcar_experiment_config(),
  backend = c("cpp", "R", "reference"), smoke_test = FALSE)
}
\arguments{
\item{p_list}{Dimensions of the data views.}
\item{n_grid}{Sample sizes, taken as nested prefixes of one simulated realization.}
\item{rank_grid}{Latent ranks.}
\item{...}{Named overrides of configuration fields. Inspect \code{names(egcar_experiment_config())}.}
\item{output_dir}{Destination for checkpoints, CSV diagnostics, and plots. Existing files of the same names are overwritten. Use separate output directories for separate runs.}
\item{workers}{Common fold-level worker budget, capped at the number of folds.}
\item{n_reps}{Number of repetitions.}
\item{config}{Configuration list returned by \code{egcar_experiment_config()}.}
\item{backend}{Compiled, accelerated R, or dense-reference EGCAR implementation.}
\item{smoke_test}{Use small grids and short iteration caps for a diagnostic run. Requires the default view dimensions and is not a performance study.}
}
\details{Exactly four EGCAR curves are produced: L11-CV, L11-rate, L21-CV,
and L21-rate. Both CV grids default to \code{10^seq(-5, 4)} and reject zero.
The L21 rate constant \code{rate_c_g} defaults to the independent pilot-calibrated
coefficient used by \code{egcar_rate(penalty="l21")}. Explicit overrides remain
effective. Invalid requested-rank loadings are recorded as
\code{status="invalid_loading"}, even when ADMM converged. Results include
\code{loading_valid}, \code{selected_rows} (for operator estimates),
\code{zero_solution}, \code{rate_zero_bound} and \code{rate_zero_certified}.
The last two fields apply to L21-rate. Unavailable subspace errors remain NA.
SGCA, RGCCA, SGCCA, MultiCCA, Oracle1-population, and Oracle2-support are retained.
SGCA initialization is bundled. The configuration field \code{sgca_time_limit}
sets one elapsed-time budget in seconds for each SGCA CV-and-refit run; its
default is 21600 (six hours). SGCA defaults to 1000 initializer iterations at tolerance 0.005 and 15000 fixed TGD updates; see comparison_cv. RGCCA and PMA are needed for the full comparison;
ggplot2 is needed for plots. Future and future.apply enable parallel CV.
The population oracle retains the original zero-penalty consensus algorithm;
neither its data nor its coefficients depend on sample size. The support
oracle retains its sample-covariance and known-support definition.
All comparison methods use the supplied local common-loss CV routines, not
alternative CV functions from their underlying libraries.
Package loading/native compilation and plot generation are outside method fit
and tuning times. All-candidate CV failure is recorded explicitly rather than
substituting a coefficient outside the supplied grid.}
\value{The configuration function returns a validated list. The runner invisibly
returns results, candidate and fold tables, fitted objects, populations,
metadata, and the path of \code{egcar_simulation.rds}.}
\examples{
cfg <- egcar_experiment_config()
stopifnot(all(cfg$rho_e_cv_grid > 0), all(cfg$lambda_g_cv_grid > 0))
\dontrun{
run_egcar_experiments("full_study", workers = 5, n_reps = 1)
}
}

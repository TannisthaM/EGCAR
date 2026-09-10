# Repository organization

`R/egcar.r`: generic EGCAR public functions.
`R/egcar_methods.R`: explicit L11/L21 fit, CV, and rate wrappers.
`R/reduced_rank_regression.R`, `R/group_reduced_rank_regression.R`: accelerated solvers.
`R/reference.R`: dense reference backends and private zero-penalty oracle.
`R/alt_SGCA.R`, `R/alt_RGCCA.R`, `R/alt_SGCCA.R`, `R/alt_MultiCCA_CrossValidation.R`: local comparison CV.
`R/sgca_initializer.R`: the bundled four-function initializer dependency closure.
`R/experiment_config.R`: shared full-study defaults and strict positive-grid validation.
`R/experiment_engine.R`, `R/experiment_loop.R`, `R/experiments.R`: full-study orchestration and plots.
`src/`: compiled native solver and registered Rcpp interface.
`inst/examples/`: complete runnable examples.
`inst/standalone/`: matching independent local script; not sourced by package functions.
`tests/testthat/`: package regression tests.
`tools/`: installation/checking/equivalence utilities.

For GitHub Desktop, select the extracted `egcar` directory, not its parent or the ZIP.
At the published repository root, DESCRIPTION and R/ must appear directly.

# Repository organization

`R/egcar.r`: generic EGCAR public functions.
`R/egcar_methods.R`: explicit L11/L21 fit, CV, and rate wrappers.
`R/reduced_rank_regression.R`, `R/group_reduced_rank_regression.R`: accelerated solvers.
`R/l21_direct.R`: current global row-L21 formulation and native-coordinate interface.
`R/reference.R`: dense reference backends and private zero-penalty oracle.
`R/alt_SGCA.R`, `R/alt_RGCCA.R`, `R/alt_SGCCA.R`, `R/alt_MultiCCA_CrossValidation.R`: local comparison CV.
`R/sgca_initializer.R`: the bundled four-function initializer dependency closure.
`R/experiment_config.R`: shared full-study defaults and strict positive-grid validation.
`R/experiment_engine.R`, `R/experiment_loop.R`, `R/experiments.R`: full-study orchestration and plots.
`src/`: compiled native solver and registered Rcpp interface.
`inst/examples/`: short `quickstart.R` and `compare_methods.R`, plus larger historical examples.
`inst/standalone/`: historical independent local script for numerical comparison; not sourced by package functions and not updated with 0.2.11 optimizations.
`tests/testthat/`: package regression tests.
`tools/`: installation/checking/equivalence utilities.

For GitHub Desktop, select the extracted `egcar` directory, not its parent or the ZIP.
At the published repository root, DESCRIPTION and R/ must appear directly.

The top-level README describes the current API. Historical algorithm and optimization reports live in `inst/doc/`; release history is in `NEWS.md`.

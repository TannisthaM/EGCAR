# egcar 0.2.3

* Memory-first experiment defaults: full fitted objects, fold-level CV tables, full loading-plot data, and per-configuration loading PDFs are no longer retained/generated unless explicitly requested.
* Save compact per-configuration loading archives with xz compression so loading heatmaps can be reconstructed without complete solver objects.
* Drop unused fold indices/training means and dispatch fold objects directly to workers to reduce future-global export pressure.
* Compact final experiment ADMM results: rate fits and selected refits keep C_hat and diagnostics but discard warm-start auxiliary state.
* For group ADMM, compute active rows directly from endpoint copies and avoid constructing view-wide G/V matrices unless state retention is requested.
* Compact SGCA caches by retaining rank-r U plus scalar initializer diagnostics rather than Pi/H/Gamma; remove an unused cached U0 matrix.
* Suppress complete external benchmark fit objects in the experiment by default; public comparison calls retain their previous behavior.
* Reduce the default loading-factor cache from 128 to 4 entries, default `keep_full_C=FALSE` because edge blocks already represent the operator, construct projector differences on demand, stop writing duplicate plot-data CSVs, and use xz-compressed checkpoints.
* Numerical estimators, tuning grids, folds, seeds, validation losses and parameter-selection rules are unchanged.

# egcar 0.2.2

* Explicitly construct each native state block as an R numeric matrix, including
  rectangular, one-row, one-column, and 1-by-1 blocks.
* Check the native matrix API version and returned state shapes before computing
  row norms or extracting loadings; do not silently reshape dimensionless arrays.
* Apply the same serialization fix to the bundled standalone runner and use a
  new native-cache filename/option there.
* Add installed matrix-interface, rank 1/2/5, warm-start, rate and CV regression
  checks, plus focused testthat coverage.
* Leave numerical solver updates, penalties, grids, seeds, benchmarks, oracles,
  experiment defaults and parallel-worker allocation unchanged.
* R/Rcpp execution of these new checks has not been performed in the patching
  environment. See inst/doc/MATRIX_INTERFACE_FIX_022.md for validation scope.

# egcar 0.2.1

* Speed-only revision of the compiled EGCAR L11/L21 backend.
* Adapted mapped Eigen products for small matrices after reviewing EfficientCCA.
* Kept Armadillo kernels for tiny and large products; disabled Eigen threading.
* Removed copies of read-only double-valued native context matrices.
* Reused projection, lifting and group-norm buffers; fused proximal/dual/residual loops.
* Added RcppEigen as a build dependency, without requiring EfficientCCA or SMUT.
* Preserved every R implementation, original controls, other methods, CV and oracles.
* Added a one-configuration, ten-method check using five shared folds/workers.
* Added native old/new comparisons, standalone timing checks and R-level tests.

# egcar 0.2.0

* Separate explicit EGCAR_L11 and EGCAR_L21 fit/CV/rate entry points.
* Strictly positive EGCAR CV grids, including the full and smoke experiments.
* Removed the redundant tied experiment and legacy result names.
* Bundled the four required SGCA initializer functions and retained MIT notices.
* Retained the local experiment's SGCA/RGCCA/SGCCA/MultiCCA CV routines.
* Added a shared-config full experiment runner and matching standalone source.
* Added explicit CV-failure diagnostics, positive-grid and initializer tests.
* Fixed the stats::toeplitz namespace import missing in the prior package.
* Retained the accelerated C++ core, dense references, and original oracles.

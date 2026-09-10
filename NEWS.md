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

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

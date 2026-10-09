# Examples

Start with **[quickstart.R](quickstart.R)**: a small simulation, six EGCAR
fits, a status/timing table, and loading extraction. It needs only the
package's required dependencies and defaults to one worker.

```r
source(system.file("examples", "quickstart.R", package = "egcar"))
```

**[compare_methods.R](compare_methods.R)** shows the four rate/CV EGCAR
variants, SGCA, RGCCA, SGCCA, MultiCCA, and two truth-informed oracles.
Install `RGCCA` and `PMA` first. Every method uses the same simulated dataset;
CV wrappers reuse the same folds. The demonstration uses small tuning grids
and a 120-second whole-CV SGCA budget, with numerical defaults preserved.
Timeouts and nonconvergence remain visible in its report.

```r
source(system.file("examples", "compare_methods.R", package = "egcar"))
```

Both scripts can also be run from a cloned repository, for example
`Rscript inst/examples/quickstart.R`, after installing `egcar`.

## Additional examples and checks

| File | Purpose |
|---|---|
| [00_install_dependencies.R](00_install_dependencies.R) | Explicit dependency installation |
| [01_quickstart.R](01_quickstart.R) | Original four-variant EGCAR example |
| [02_all_methods_cv.R](02_all_methods_cv.R) | Full original local study; **not a short demonstration** |
| [03_parallel_cv.R](03_parallel_cv.R) | Fold-level parallel CV |
| [04_backend_validation.R](04_backend_validation.R) | Backend comparison |
| [05_backend_timing.R](05_backend_timing.R) | Timing on your R installation |
| [06_comparison_cv.R](06_comparison_cv.R) | Four comparison wrappers; SGCA uses its current score default |
| [06_small_all_methods_check.R](06_small_all_methods_check.R) | Ten-method configuration with five folds/workers |
| [07_matrix_interface_check.R](07_matrix_interface_check.R) | Matrix input checks |
| [08_performance_check.R](08_performance_check.R) | Larger numerical/performance diagnostics |

The historical numbered scripts retain their original study settings.
Use the short examples above for first use. Fitting functions never install
optional packages automatically. `workers` distributes folds within a method;
it does not submit methods as separate jobs.

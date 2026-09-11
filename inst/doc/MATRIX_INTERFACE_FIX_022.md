# Matrix-interface patch: egcar 0.2.2

## Basis and diagnosis limits

This patch is based on the uploaded `egcar_0.2.1_github(1).zip`, not an independently
verified newer GitHub checkout. The reported error is:

```
Error in base::rowSums(x, na.rm = na.rm, dims = dims, ...) :
  'x' must be an array of at least two dimensions
```

In the supplied source, `src/egcar_native.cpp` returns the solver's nested
`std::vector<arma::mat>` state fields through generic `Rcpp::List::create`
conversion. The accelerated R dispatcher does not check the returned matrix
shapes before the L11/L21 wrappers and loading routines use them. In particular,
`egcar_loading_from_operator()` calls `egcar_group_norms()`, which calls
`rowSums()` and `colSums()` on individual edge blocks.

That return boundary is the suspected source of dimension loss. No traceback
from the original package failure was supplied, and R/Rcpp are not installed
in the patching environment. Therefore this is a targeted interface fix, not
a claim that the original failure has been reproduced or that every version
of Rcpp necessarily flattens such containers. The included installed-package
check tests the actual boundary and numerical results on the user's system.

## Implementation

| File | Change |
|---|---|
| `src/egcar_native.cpp` | Allocate an explicit `Rcpp::NumericMatrix` for every returned C/Z/H or C/Gk/Gl/Vk/Vl block; copy the original column-major values; return an internal `matrix_api = 2L` marker. |
| `R/utils.R` | Add `.egcar_check_solver_state()`, validating the required state lists and each expected `p_k x p_l` matrix shape. |
| `R/reduced_rank_regression.R` | Check the compiled API marker and state dimensions immediately after the solver returns, before sparsity postprocessing and loading extraction. |
| `inst/standalone/run_EGCAR_local.R` | Apply the same explicit conversion and checks to the separately shipped standalone source; change its native cache filename and runtime-option key. |
| `inst/validation/matrix_interface.R` | Add direct native-interface checks and public fit/path/rate/CV regression checks. |
| `inst/examples/07_matrix_interface_check.R` | Run the installed check with five CV workers by default and save diagnostic CSV/RDS files. |
| `tests/testthat/test-native-interface-022.R` | Test malformed-state rejection, stale-DLL rejection, and the installed serial regression suite. |
| `DESCRIPTION`, `NEWS.md`, `README.md` | Set the local patch version to 0.2.2 and document the changes and validation limits. |

The new native helper returns owning R matrices with explicit dimensions,
including 1-by-1, 1-by-p, p-by-1 and rectangular blocks. No transposition is
introduced. The R checker does **not** call `as.matrix()` on dimensionless
vectors or guess a missing shape: it rejects malformed output with the field,
edge index, expected dimensions, and observed type/shape.

The raw API marker identifies mismatched R code and an older loaded DLL.
Updating source files without reinstalling and restarting R is not sufficient
for a compiled-package change. The registered native function signature is
unchanged, so neither `R/RcppExports.R` nor `src/RcppExports.cpp` is modified.
The marker is internal; the public fit/CV return structure is not redesigned.

## What is unchanged

The native numerical core is byte-for-byte unchanged. The dense-reference and
optimized-R arithmetic, ADMM objective, proximal rules, residual calculations,
adaptation rules, tolerance/iteration defaults, CV grids, warm-start convention,
seeds, covariance preparation, loading formulas, metrics, benchmark routines,
and oracles are retained.

`R/experiment_config.R` is unchanged. Its existing defaults already use
signal 0.8, ranks 1/2/5, three views of dimension 15, and sample sizes
30, 45, 100, 1000, 5000 and 10000. Worker allocation is not changed: specify
`workers = 5L` when running the package experiment. Validation of matrix
shapes adds a small once-per-solve check, not another ADMM iteration or penalty.

The separately bundled standalone runner retains its own previous numerical
implementation. This patch does not silently replace it with a different
version of the accelerated core.

## Validation actually executed in the patching environment

System Armadillo and Eigen headers were used to compile and run the repository's
standalone native tests with a C++14 compiler. These tests exercise the native
arithmetic, **not Rcpp conversion, package installation, R code, or CV**.

* Dense-reference comparison: 120 cases passed. Maximum row-sum norm of the
  C-state difference was approximately `8.83e-11`. Iteration counts, convergence
  flags and augmentation parameters matched.
* Comparison with the previous native core: 864 cases passed. Maximum absolute
  state difference was approximately `1.07e-14`; stopping/adaptation decisions
  matched.
* Native mapped-context checks: 8 cases passed. Read-only inputs were preserved.

The recorded results and scope are in
`inst/validation/matrix_patch_022_native_checks.json`. A source comparison also
confirmed that the native core and the unchanged implementation files were not
modified. The patch was checked for clean application to the uploaded source.

**Not executed here:** R parsing via `parse()`, source-package installation,
Rcpp-wrapper compilation, the installed regression check, `testthat`,
`R CMD check`, five-worker CV, the full all-method experiment, or macOS testing.
No statement above should be interpreted as an R-level test pass.

## Checks to run after installation

Restart R before the clean source installation, and use the complete package
archive rather than treating the changed-files-only archive as a package.
After installation, `packageVersion("egcar")` should report `0.2.2`.

```r
source(system.file("examples", "07_matrix_interface_check.R",
                   package = "egcar", mustWork = TRUE))
```

The default check uses five CV workers and requires `future` and `future.apply`.
It does not need RGCCA, PMA, ggplot2, or testthat. It includes 16 direct
native case/family combinations covering one-row/one-column blocks, rectangles,
wide data, an empty spectral basis, and the user's dimensions/ranks. It checks
raw shapes before the public wrappers operate on the results, compares native
state values with the optimized-R backend, and checks input preservation.

At signal 0.8 and ranks 1, 2 and 5, it also checks both penalty families' fixed
fits, warm starts, regularization paths, rate fits, and serial versus
five-worker CV on the same five folds and two small positive candidates.
The rate penalty is allowed to produce a legitimate `invalid_loading` result
when it selects too few rows; that is reported, not turned into a fabricated
rank-r loading. Nonconvergence under diagnostic iteration limits is likewise
reported separately from matrix-interface success.

For a serial-only check:

```r
Sys.setenv(EGCAR_MATRIX_CHECK_WORKERS = "1")
source(system.file("examples", "07_matrix_interface_check.R",
                   package = "egcar", mustWork = TRUE))
Sys.unsetenv("EGCAR_MATRIX_CHECK_WORKERS")
```

The normal five-worker example saves reports under
`egcar_matrix_interface_check/`. An alternate location can be selected with
`EGCAR_MATRIX_CHECK_OUTPUT`. The check restores a pre-existing future plan.

The original all-ten-method example is unchanged and remains available:

```r
source(system.file("examples", "06_small_all_methods_check.R",
                   package = "egcar", mustWork = TRUE))
```

For the complete package check, run `Rscript tools/check.R .` from the repository
root after installing the check dependencies. The existing GitHub Actions
workflow will run the new testthat file when the patch is committed and pushed;
this patch delivery itself does not modify the remote repository.

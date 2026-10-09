# EGCAR

**Sparse generalized correlation analysis for multiple data views, through regression.**

EGCAR finds shared low-dimensional structure in datasets measured on the same
observations, such as several molecular assays on the same samples. It estimates
a sparse cross-view regression operator, then extracts loading vectors that map
each view to shared components. The R package provides separate entrywise
(`L11`) and row-group (`L21`) penalties, cross-validation, rate-scaled penalties,
and wrappers for four comparison methods.

[Quick start](#quick-start) · [Methods](#methods) ·
[Examples](inst/examples/README.md) · [Release notes](NEWS.md) ·
[Issues](https://github.com/TannisthaM/EGCAR/issues)

## Installation

Install from GitHub in a fresh R session:

```r
install.packages("remotes")  # Once, if needed.
remotes::install_github("TannisthaM/EGCAR", dependencies = NA,
                        upgrade = "never", build_vignettes = FALSE)
library(egcar)
packageVersion("egcar")
```

R 4.0 or later and a C++14 toolchain are required. Source installation needs
Rtools on Windows, the R-compatible compiler toolchain on macOS, or the R
development/build tools on Linux. Required R dependencies are installed by
`remotes`; optional comparison and parallel packages can be added separately:

```r
install.packages(c("RGCCA", "PMA"))            # Optional comparison methods.
install.packages(c("future", "future.apply")) # Optional parallel CV folds.
```

The optional packages' own R-version requirements also apply.

For a downloaded source tarball with dependencies already installed:

```r
install.packages("egcar_0.2.18.tar.gz", repos = NULL, type = "source")
```

The repository itself is the package root: `DESCRIPTION`, `R/`, and `src/`
belong directly under `EGCAR/`. There is no additional nested package folder.

## Quick start

Run four EGCAR variants on one small simulated dataset. Preparing the folds
once ensures that both CV methods use exactly the same train/validation splits.

```r
library(egcar)

sim <- egcar_simulate(
  n = 90, p_list = c(6, 6, 6), rank = 1, active_per_view = 3,
  toeplitz_rho = c(0.5, 0.7, 0.9), signal = 0.8, seed = 101
)
shared <- egcar_cv_data(sim$views, nfolds = 3, seed = 102)
grid <- c(0.001, 0.01, 0.1)  # Small example grid; choose for your application.

fits <- list(
  L11_rate = EGCAR_L11_Rate(shared, rank = 1),
  L21_rate = EGCAR_L21_Rate(shared, rank = 1),
  L11_CV = EGCAR_L11_CV(shared, rank = 1, lambda = grid),
  L21_CV = EGCAR_L21_CV(shared, rank = 1, lambda = grid)
)
invisible(lapply(fits, print))

fit <- fits$L21_CV
if (identical(fit$status, "ok") && isTRUE(fit$converged)) {
  loadings <- coef(fit, by_view = TRUE)  # One p_k x rank matrix per view.
  scores <- predict(fit, newdata = sim$views, type = "views")
  print(lapply(loadings, dim))
}
```

To choose a penalty coefficient directly, use
`EGCAR_L11(shared, rank = 1, lambda = 0.01)` or
`EGCAR_L21(shared, rank = 1, lambda = 0.01)`.
The generic equivalents are `egcar_fit()`, `egcar_rate()`, and `egcar_cv()`
with `penalty = "l11"` or `penalty = "l21"`.

The [short runnable script](inst/examples/quickstart.R) includes all six calls,
a status/timing table, and loading extraction. After installation:

```r
source(system.file("examples", "quickstart.R", package = "egcar"))
```

## Your own data

Supply a list of at least two finite numeric matrices. Rows must refer to
the **same observations in the same order**; columns are view-specific features.
The number of columns may differ between views.

```r
# X1, X2 and X3 are your matched sample-by-feature matrices.
views <- list(view1 = X1, view2 = X2, view3 = X3)
shared <- egcar_cv_data(views, nfolds = 5, seed = 2026)
fit <- EGCAR_L21_CV(shared, rank = 2, lambda = c(0.001, 0.01, 0.1))
```

The package centers each training split and applies its training means to
validation data. It does not automatically standardize variables or blocks,
impute missing values, or align sample IDs. Choose any scaling deliberately;
for evaluation, estimate preprocessing parameters using training data only.
Prediction requires all views and the original feature names/order.

## Methods

| Method | Public function | Penalty or tuning |
|---|---|---|
| EGCAR L11 | `EGCAR_L11()` | User-specified entrywise penalty |
| EGCAR L21 | `EGCAR_L21()` | User-specified global row-group penalty |
| EGCAR L11, rate | `EGCAR_L11_Rate()` | Rate-scaled coefficient |
| EGCAR L21, rate | `EGCAR_L21_Rate()` | Rate-scaled coefficient |
| EGCAR L11, CV | `EGCAR_L11_CV()` | Direct coefficient grid |
| EGCAR L21, CV | `EGCAR_L21_CV()` | Direct coefficient grid |
| Sparse GCA (Gao and Ma) | `sgca_cv()` | Total row sparsity by default |
| RGCCA | `rgcca_cv()` | Regularization parameter; needs `RGCCA` |
| SGCCA | `sgcca_cv()` | Sparsity bound; needs `RGCCA` |
| MultiCCA | `multicca_cv()` | L1 bound; needs `PMA` |
| Population oracle | `egcar_oracle_population()` | Known population covariances |
| Support oracle | `egcar_oracle_support()` | Known true variable supports |

The oracles are simulation diagnostics: their population information is not
available in an ordinary data analysis. The
[ten-method example](inst/examples/compare_methods.R) combines the four
rate/CV EGCAR variants, four comparison methods, and two oracles on the same
dataset and folds. It uses small grids and a 120-second whole-CV SGCA budget
for demonstration; package defaults remain unchanged.

```r
source(system.file("examples", "compare_methods.R", package = "egcar"))
```

SGCA uses its paper covariance-trace validation score by default; the other
CV methods use the common generalized Rayleigh score. Explicitly setting
`benchmark_control(sgca_cv_score = "common_loss")` selects the common-score
SGCA variant. Shared folds alone do not establish equal numerical accuracy.
See the [SGCA specification](inst/doc/SGCA_PAPER_0216.md) for defaults and diagnostics.

## How the method works

EGCAR fits cross-view coefficient blocks using a convex quadratic regression
loss and **one** sparsity penalty at a time. L11 encourages individual zero
coefficients. The current L21 penalty is

$$
  \lambda\lVert C\rVert_{2,1}
  = \lambda\sum_i\left(\sum_j C_{ij}^2\right)^{1/2}.
$$

In version 0.2.18, L21 is optimized over directed off-diagonal blocks using
ADMM. Its sparse estimate is then symmetrized as `(C + t(C)) / 2`, before
extracting the requested loading vectors. Diagonal blocks remain zero.
Convergence describes the directed optimization; this final projection is
not claimed to solve the earlier symmetry-constrained objective and can
expand row support. The inherited L21 rate multiplier, `0.05`, has not been
recalibrated for this formulation.

The implementation reduces repeated work through cached covariance/spectral
quantities, compact native L21 states, warm starts along CV paths, and partial
eigenvalue decomposition for loading extraction. CV can distribute folds
across workers. Runtime and peak memory depend on the data and iteration count;
there is no universal speedup guarantee. See the
[implementation notes](inst/doc/L21_IMPLEMENTATION_0218.md) and
[measured validation](inst/doc/VALIDATION_0218.md).

## Parallel cross-validation

```r
# Requires future and future.apply; use a prepared three-fold object.
fit_parallel <- EGCAR_L21_CV(shared, rank = 1,
                            lambda = c(0.001, 0.01, 0.1), workers = 3)
```

`workers` parallelizes **folds within one method**, up to the number of folds.
Candidates within each EGCAR fold run sequentially with warm starts. Calling
several methods in a list still runs those methods sequentially; separate
Slurm jobs are a distinct experiment-launching choice. The default is one
worker. Match workers to the CPUs allocated to your R process.

## Results, convergence, and timing

Inspect `fit$status`, `fit$converged`, and, for CV, `fit$cv_fold_table` before
using results. `coef(fit)` returns the stacked loading matrix, while
`coef(fit, by_view = TRUE)` splits it by view. `predict()` applies the stored
training means; `egcar_score(fit, newdata)` evaluates the generalized
Rayleigh score. The quick-start prediction above is in-sample, not an estimate
of test-set performance.

| Field | Meaning |
|---|---|
| `fit_time` | Complete final fit, through calculation of the loading vectors |
| `solver_time`, `loading_time` | EGCAR components already included in `fit_time` |
| `tuning_time` | CV fitting and scoring, including worker overhead when applicable |
| `total_time` | For EGCAR CV, `tuning_time + fit_time` |
| `wall_time` | Current public call, including preparation performed in that call |
| `preparation_time` | Preparation cost associated with the input; it may have been incurred earlier |

Do not add `loading_time` to `fit_time` again. When reusing `shared`, account
for its preparation once. SGCA's finite-budget supervisor also tracks elapsed
time across its subprocess; see its method-specific documentation.

An iteration cap, a timeout, and convergence are different outcomes. A fit
without the requested positive loading rank is reported as `invalid_loading`;
CV with no eligible candidate is `no_valid_cv`. SGCA's fixed-iteration
completion can have `status = "ok"` while `converged = FALSE`. The package
does not silently lower the requested rank to make a failed fit usable.

## Documentation and development

```r
help(package = "egcar")
?egcar_fit
?egcar_cv
?egcar_control
?comparison_cv
```

- [Example index](inst/examples/README.md): short examples and larger studies.
- [Repository guide](inst/doc/REPOSITORY_LAYOUT.md): where the code lives.
- [Release notes](NEWS.md): version history and behavior changes.
- [Packaging validation](inst/doc/PACKAGE_LAYOUT_0218.md): checks for this documentation update.
- [Contributing](CONTRIBUTING.md): source installation, checks, and issue reports.

EGCAR is distributed under [GPL-3](LICENSE.md). Bundled initializer helpers
retain their separate MIT notice in [COPYRIGHTS](inst/COPYRIGHTS).

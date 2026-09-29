# SGCA specification in EGCAR 0.2.16

This release follows **the paper's Algorithm 1**, using Sections 5.1–5.2 for the
standard sparse-GCA benchmark. It is not a byte-for-byte port of the authors'
R script, which differs from the printed algorithm in several places. It does
not reproduce every simulation design in the paper or establish matched accuracy
between SGCA and EGCAR.

## Sources and scope

- Gao and Ma (2023), *Sparse GCA and Thresholded Gradient Descent*, JMLR 24(135):
  [paper](https://www.jmlr.org/papers/volume24/21-0745/21-0745.pdf).
- [Authors' code](https://github.com/sggao/sparse-gca/tree/a121e7f3bb72f825117c11a6d0bb3273bcb2bc81),
  inspected at commit `a121e7f3bb72f825117c11a6d0bb3273bcb2bc81`.
- Initializer dependency originally bundled from
  [TannisthaM/SGCA](https://github.com/TannisthaM/SGCA/blob/main/R/gao_cv_functions.R),
  with its MIT notice retained. The helper is now adapted to the paper constraint.

The paper minimizes `-tr(V' A V) + lambda/2 * ||V' B V - I||_F^2`, where
`A` is the joint covariance and `B` its block diagonal. Its initialization uses
trace-equality Fantope projection, leading eigenvectors and row truncation.
Algorithm 1 normalizes and rescales the start, repeatedly applies a gradient step
and row truncation, then normalizes the final loading. Section 5.1 uses sparsity
20; Section 5.2 selects sparsity by five-fold CV over 5,10,...,100.

## Effective defaults and their provenance

| Setting | Package default | Provenance |
|---|---|---|
| Covariance divisor | `n_train` | Paper Algorithm 1 discussion |
| Added metric ridge | `0` | No ridge in Algorithm 1 |
| Initializer penalty | `0.5*sqrt(log(p_total)/n_train)` | Sections 5.1–5.2 |
| Fantope constraints | eigenvalues in [0,1], trace=rank | Equation (14) |
| Initializer change tolerance | `0.005`, absolute Frobenius matrix change | Authors' R default |
| Initializer outer cap | `1000`, actually enforced | Authors' advertised R default, counter bug repaired |
| Initializer inner updates | `20`; ADMM nu=1 | Authors' R implementation |
| TGD lambda | `0.01` | Section 5.1 benchmark |
| TGD eta | `0.001`; update uses `2*eta` times gradient-half | Algorithm 1 and Section 5.1 |
| TGD updates | exactly `15000` if successful | Sections 5.1–5.2 |
| TGD change diagnostic | `1e-6`, absolute Frobenius change | Authors' R example; diagnostic only in fixed-step mode |
| Sparsity grid | valid members of 5,10,...,100 | Section 5.2 |
| CV score | test covariance trace, maximized | Section 5.2 |
| Whole CV+refit time budget | six hours | EGCAR operational policy, not a paper prescription |

The initializer tolerance is an implementation stopping test, not a guarantee
that the convex program is solved to a particular statistical error. The paper
does not prescribe this numeric tolerance. The reference code's MATLAB
initializer has different controls and stopping criteria; they are not combined
with the R tolerance in this release.

## Changes from EGCAR 0.2.15

1. Replace the scaled-iterate/backtracking routine with the unscaled paper update.
   Initialization is normalized in the training metric before rescaling.
   No line search or silent step-size reduction is performed.
2. Replace the initializer's trace-at-most-r projection with trace-equals-r.
   Projection uses an eigenvalue clipping shift found by bounded bisection.
   Leading algebraic eigenpairs replace SVD for initializer loading extraction.
3. Apply the rate-scaled initializer penalty independently on training folds and
   the full-sample refit. Here n means the number of observations used in that
   fit; this is our fold-wise application of the paper rate, not a claim about
   unpublished CV implementation details. Cache one initializer per fold and actual rho value.
   Numeric `rho_grid` still means direct coefficients for backward compatibility.
4. Tune sparsity alone by default. Explicit rho/lambda grids remain supported,
   but are custom experiments. Sparsity always counts rows across all views.
5. Preserve the EGCAR solvers, their convergence defaults, and the complete fit
   timing through loading extraction introduced in 0.2.15.

The authors' R TGD uses eta rather than 2*eta and omits Algorithm 1's initial
metric normalization; it also uses `cov()` in the example and SVD extraction.
Its rank-greater-than-one threshold helper indexes the matrix incorrectly and
relies on an external `r`. Those behaviors are not reproduced. The paper's
smaller-index tie rule is enforced. The R initializer reuses its outer counter
inside its inner loop; our separate counters enforce the cap. No external
function bodies from the authors' unlicensed repository were added to EGCAR.

## Completion and convergence

`status="ok"` and `completed=TRUE` mean the requested procedure completed and
returned finite, full-rank loadings. Under fixed-step benchmarking, this does
**not** imply `converged=TRUE`.

A finite initializer that reaches its cap can proceed to TGD. This makes the
configured finite initializer budget an explicit approximation. Its
`sgca_init_converged` flag remains FALSE and the cap is recorded. Set
`sgca_cv_require_convergence=TRUE` to require both numerical change criteria on
all training folds and the refit. No candidate with failed loading extraction,
non-finite iterates or an incomplete TGD run is eligible.

A finite time limit still covers SGCA preparation, the entire CV grid and final
refit. Shared data/fold preparation occurs before that budget. On timeout, no
completed estimate is returned and the status remains `time_limit`. The SLURM
job wall time is separate and must allow the complete requested cluster task.

Diagnostics are retained even when dense benchmark fits are dropped:

- `fit$diagnostics`: actual full-sample rho; initializer and TGD iteration counts,
  change magnitudes, convergence flags and stop reasons; algorithm version.
- `fit$cv_fold_table`: the same fields for individual candidates and folds,
  together with finite-score and completion flags.
- Study `simulation_results.csv` and, when enabled, `external_cv_fold_results.csv`
  carry these diagnostics. `benchmark_completed` is independent of `converged`.
- `fit_time` includes final loading normalization; `tuning_time` includes fold
  fits and SGCA preparation/caching. `total_time=fit_time+tuning_time`.

## Use the defaults

```r
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.16")
shared <- egcar_cv_data(views, nfolds=5, seed=1)
fit <- sgca_cv(shared, rank=1)
fit$status
fit$converged
fit$diagnostics
fit$cv_table
fit$cv_fold_table
```

For the Section 5.1 fixed sparsity setting, supply `k_grid=20` (requires p>=20).
That single candidate is still evaluated on the supplied folds before refitting.
For very small problems with no valid default grid member, the public default
uses k=p; this is a dimension accommodation, not a paper benchmark setting.

Optional common-score comparison, keeping the paper solver:

```r
fit_common <- sgca_cv(shared, rank=1,
  benchmarks=benchmark_control(sgca_cv_score="common_loss"))
```

Optional numerical-tolerance experiment:

```r
fit_tol <- sgca_cv(shared, rank=1,
  benchmarks=benchmark_control(sgca_stopping="absolute_change",
    sgca_tgd_tol=1e-6, sgca_tgd_max_iter=15000,
    sgca_cv_require_convergence=TRUE))
```

The second variant permits early stopping and is not the fixed-step benchmark.
With `sgca_stopping="absolute_change"`, `sgca_tgd_max_iter=Inf` is permitted; the
separate elapsed budget still applies. Fixed-iteration mode rejects an infinite
iteration count. `egcar_score()` always reports the common score, which differs
from SGCA's default selection score. Do not compare those CV-score columns as
though they measured the same quantity.

## Memory and runtime

Initializer matrix factors and scaled constants are prepared once per training
sample. The initializer is reused across sparsity candidates, and normalized
starts can be reused across lambda values. Hard thresholding and TGD operate on
p-by-r matrices. Validation contracts the held-out observations to r score
columns, avoiding a new dense covariance when the existing fold representation
stores observations. Dense fold covariances now retain their p-length mean
vector so test-centered scoring can be recovered without retaining raw views.

The initializer still stores dense p-by-p matrices and performs full
eigendecompositions. This is a substantial high-dimensional cost. Fixed-step
TGD performs all requested updates even when the diagnostic change becomes
small. This release does not claim a general runtime speedup or cluster-scale
memory measurement. Equal iteration counts or tolerance numbers do not establish
matched statistical accuracy with EGCAR.

## Upgrade

Install the source package, restart R and all workers, and recreate saved CV-data
objects. Old dense validation objects do not contain the means needed for paper
scoring; the new implementation requests that they be recreated.

Copy the 0.2.16 cluster launchers into your experiment directory and choose a
fresh output directory. Existing launchers contain explicit old SGCA grids,
ridge and iteration caps that override package defaults. Old saved configuration
lists likewise retain explicit fields and are not silently rewritten. Prefer
`egcar_experiment_config()` to create a fresh configuration. The historical
standalone script now delegates to the installed package instead of maintaining
a second solver implementation.

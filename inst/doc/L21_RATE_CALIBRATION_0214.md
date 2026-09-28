# EGCAR-L21-rate: calibration and correction in 0.2.14

## Result and scope

Version 0.2.14 changes the default L21 rate multiplier from 1 to **0.05**.
It also corrects experiment status reporting: convergence of the coefficient
optimizer no longer causes an unavailable loading to be labeled a success.
These are a penalty calibration and a diagnostic correction, not a change to
the L21 objective or a proof of successful high-dimensional sparse recovery.

The selected multiplier produced valid, converged requested-rank loadings in
all 36 pilot cases; the old multiplier did so in 4. However, the selected fits
retained every variable. The choice optimizes the prespecified held-out score
over a finite grid, not support recovery. The raw results and reproducible
scripts are included so that this distinction can be assessed directly.

## Why zero estimates can be legitimate optima

For centered training views, let `S_kl = crossprod(X_k, X_l) / n`. The package
uses the unchanged penalty rule

\[
\lambda_g=c_g\sqrt{\frac{d_{\max}+\log p}{n}},\qquad
d_{\max}=\max_k(p-p_k),\qquad p=\sum_k p_k.
\]

For three views of 1,000 variables and `n=150`, the old rule gives
`lambda_g = 3.658785` and the new default gives `lambda_g = 0.182939`.
The latter is twenty times smaller. The covariance divisor, overlapping group
definition and objective scaling have not been changed to disguise this choice.

Write `M_k(C)` for the block row containing all edges incident to view `k`,
transposing an edge when necessary. Each edge entry belongs to both endpoint
groups. The minimized objective is

\[
F(C)=\sum_{k<l}\left\{\tfrac12
\operatorname{tr}(C_{kl}^{T}S_{kk}C_{kl}S_{ll})-
\langle S_{kl},C_{kl}\rangle\right\}
+\lambda_g\Omega(C),\qquad
\Omega(C)=\sum_k\|M_k(C)\|_{2,1}.
\]

Define the sufficient bound

\[
b=\tfrac12\max_{k,i}\|M_k(S)_{i,:}\|_2.
\]

Counting the two endpoints and applying Cauchy--Schwarz gives

\[
\sum_{k<l}\langle S_{kl},C_{kl}\rangle
=\tfrac12\sum_k\langle M_k(S),M_k(C)\rangle
\le b\Omega(C).
\]

The quadratic terms are nonnegative for positive semidefinite within-view
sample covariances, including singular ones. Therefore

\[
F(C)-F(0)\ge(\lambda_g-b)\Omega(C).
\]

Consequently, **`lambda_g > b` certifies that zero is the unique optimum**.
The converse is not asserted: a penalty at or below the bound can still yield
an invalid or zero estimate. This is a sufficient bound, not the exact
overlapping-group critical penalty.

ADMM can correctly converge to zero. There is then no loading subspace of the
requested rank. Reporting `converged=TRUE` alongside `status="invalid_loading"`
is intentional; its subspace error remains `NA`. No arbitrary eigenvectors,
reduced rank or substituted fit are used to fill that value.

## Independent calibration

The protocol in `inst/validation/l21_rate_0214/PROTOCOL.md` was fixed before
candidate evaluation. The uploaded experiment results and their ground-truth
errors were not used to select the coefficient.

The pilot crosses four dimensions per view (10, 100, 500, 1,000), three ranks
(1, 3, 5), and three signals (0.3, 0.5, 0.8): 36 independent populations in
total. There are three views, five active variables per view, and Toeplitz
correlations 0.5, 0.7 and 0.9, following the supplied simulation model. Each
case uses 150 training and 600 independent validation observations.

Eight multipliers were tested, with a separate cold-start fit for every
non-certified candidate. Control settings were `max_iter=2000`,
`abs_tol=1e-5`, `rel_tol=1e-4`, and `check_every=5`. Certified zero candidates
were recorded as invalid without ADMM and remained in all denominators.
There were 288 pilot candidate/case combinations: 241 ADMM fits and 47
certified-zero shortcuts. The later confirmation added 75 actual ADMM fits
and 45 shortcuts. Thus 316 numerical fits were run across the two stages.

Only a multiplier with a valid, converged requested-rank loading in **every**
pilot case was eligible. Among eligible candidates, the rule maximized the
mean validation score divided by rank, with equal weight per case and exact
ties favoring the larger coefficient. This is the package's held-out
generalized Rayleigh score; it is not a probability or a support metric.
Population subspace errors and true support did not enter the selection.

| Multiplier | Valid and converged / 36 | Eligible mean score / rank | Mean selected fraction, all cases |
|---:|---:|---:|---:|
| **0.05** | **36** | **1.327040** | **1.000000** |
| 0.10 | 36 | 1.324862 | 1.000000 |
| 0.25 | 36 | 1.294710 | 1.000000 |
| 0.40 | 36 | 1.263559 | 0.998065 |
| 0.50 | 36 | 1.254431 | 0.713157 |
| 0.60 | 14 | ineligible | 0.244241 |
| 0.75 | 9 | ineligible | 0.191574 |
| 1.00 | 4 | ineligible | 0.073148 |

`0.05` is the lower boundary of the tested grid. Its small average advantage
over `0.10` does not establish a statistically significant difference or an
optimum outside that grid. The value was frozen in `selected_multiplier.txt`
before the confirmation simulations; it was not readjusted afterward.

## Independent confirmation

Confirmation used 60 **new** populations: the same dimensions and signals,
now crossed with all ranks 1 through 5. Only the frozen coefficient 0.05 and
the old coefficient 1 were evaluated. No confirmation result changed the
selection. Seeds begin in a separate range from the pilot.

| Variables per view | Cases | Valid/converged, new 0.05 | Valid/converged, old 1 | New mean score / rank | New mean subspace error / sqrt(rank) |
|---:|---:|---:|---:|---:|---:|
| 10 | 15 | 15 | 9 | 1.944374 | 0.344967 |
| 100 | 15 | 15 | 0 | 1.282764 | 0.881487 |
| 500 | 15 | 15 | 0 | 1.018171 | 0.989384 |
| 1,000 | 15 | 15 | 0 | 1.003341 | 0.997282 |
| **All** | **60** | **60** | **9** | **1.312163** | **0.803280** |

Subspace error uses the population Sigma0 metric; its normalized range is
0 to 1, with 0 representing matching subspaces. Thus the high-dimensional
values near 1 show **poor recovery despite valid loadings**. Every selected
new fit retained all variables: 30, 300, 1,500 or 3,000 rows depending on the
dimension. This calibration does not solve high-dimensional support recovery.

For the old coefficient, 45 cases were invalid by the strict zero certificate
and were not passed to ADMM. The other 15 were actually fitted; 9 were valid
and 6 invalid. All 60 new-coefficient cases were actually fitted by ADMM and
converged. Comparing a score averaged over the old coefficient's 9 survivors
with one averaged over all 60 new fits would introduce selection bias; the
summary tables label these conditional averages explicitly.

There is only one new population at each dimension/rank/signal combination
per stage. These results establish what happened in this finite study, not
a precise population failure rate or a guarantee on the original experiment
realizations. In particular, the exact uploaded tier-3 runs were not replayed.

## Implementation and memory cost

- `R/l21_rate.R` defines the shared frozen constant and the zero certificate.
  The default public API and experiment configuration both read this constant.
  `multiplier=NULL` resolves to 0.05 for L21 and 1 for L11. Explicit numeric
  values retain their previous meaning, including an explicit value of 1.
- `R/egcar.r` and `R/egcar_methods.R` expose the new default. Rate fits include
  the actual multiplier, its source, the sufficient zero bound, a certificate
  flag, and time spent computing the diagnostic.
- `R/experiment_loop.R` computes the bound for L21-rate, includes its cost in
  fitting time, and catches loading-extraction errors for that method.
- `R/experiment_engine.R` separates loading validity from optimizer convergence
  for the shared experiment evaluator. It adds `loading_valid`, `selected_rows`,
  `zero_solution`, `rate_zero_bound`, and `rate_zero_certified` to result rows.
  It preserves genuine error, timeout and skip statuses. Invalid fits enter
  the existing failure CSV, with a reason and unavailable metrics left as NA.

The certificate visits the stored cross-covariance blocks in tiles of at most
256 by 256 entries. Squared row sums contribute to the first endpoint and
squared column sums to the second. It needs `O(p + 256^2)` extra storage and
`O(sum_{k<l} p_k p_l)` work. It does not concatenate block rows or construct
a full `p` by `p` covariance. It uses the existing covariance blocks; their
storage is not removed by this diagnostic.

The public fitting function reports the certificate without bypassing the
solver or automatically reducing the penalty. Only the calibration harness
uses the proved certificate to avoid redundant fits during its grid study.
The compressed ADMM state and partial eigenpair extraction from 0.2.12 are
retained. The full experiment runner still has other dense population and
metric matrices. This release is not a claim to remove all quadratic storage.

Lowering the penalty retains more variables. It can increase loading-extraction
work and memory compared with a zero or sparse estimate. No runtime or peak-RSS
speedup is claimed for 0.2.14's calibration itself.

## Reproduction and validation

From the extracted package root, using an installed `egcar`, run:

```sh
Rscript tools/calibrate_l21_rate_0214.R pilot calibration_output
Rscript tools/calibrate_l21_rate_0214.R select calibration_output
Rscript tools/calibrate_l21_rate_0214.R confirmation calibration_output
Rscript tools/summarize_l21_calibration_0214.R calibration_output
```

The helper constructs the same latent distribution as `egcar_simulate()` but
omits an unnecessary full population eigendecomposition. A same-seed small
case was compared with the package simulator and matched to numerical
tolerance. The pilot used installed 0.2.13 with explicit candidate coefficients;
confirmation used installed 0.2.14. Their objective and numerical solvers are
identical. Seed tables and session information are shipped with the results.
For replay with a specific installation, set `EGCAR_CALIBRATION_LIBRARY` to
its R library directory. A third argument can select comma-separated scenario
indices for partitioned jobs. Use a fresh output directory for a new study;
completed per-case CSV files are treated as checkpoints.

Validation used R 4.3.3 and a compatible single-thread LP64 OpenBLAS. The BLAS
was checked against the reference implementation on matrix products and an
L21 fit. Timings in the raw CSVs include loading/metric work, omit population
and data preparation, and include certificate shortcuts. They are not an
end-to-end performance comparison or a cluster memory benchmark.

The new regression tests cover the calibrated default and explicit overrides,
the unchanged L11 default, an analytic two-scalar zero/nonzero transition,
tile boundaries and both endpoints, status propagation, and the failure CSV
from a small complete experiment. The full package check record is included
in `inst/validation/l21_rate_0214/package_check.txt`; the test output is in
`package_tests.txt`. Optional external comparison tests require their packages
and an explicit opt-in.

The completed direct test suite passed **860 expectations**, with **zero
failures, zero errors and zero warnings**; one optional comparison test was
skipped. `R CMD check --no-manual` completed with **zero errors, zero warnings
and four notes** (unavailable optional packages, the C++14 specification,
installed binary size, and the included LaTeX source). Package examples ran.
All three SLURM scripts passed `bash -n`, and their R launchers parsed.

One earlier check attempt failed at staged installation with a lazy-load
database error. A clean installation/check using the reference BLAS passed
without a package code change; the underlying cause of that earlier error was
not isolated. The check's test transcript was incomplete, so the entire suite
was also run directly against the checked installation and its final completion
marker and counts were verified. Calibration used the separate successful
OpenBLAS installation described above.

The supplied full tier-3 experiment, its original data realizations, all
comparison methods, and a six-hour SGCA run were **not** rerun for this release.
The executed pilot and independent confirmation are new simulations. This
release does not claim to fix separate comparator worker failures or out-of-
memory errors. Existing six-hour SGCA behavior from 0.2.13 is retained.

## Using the release in experiments

```r
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.14")
cfg <- egcar_experiment_config()
stopifnot(cfg$rate_c_g == 0.05)
fit <- EGCAR_L21_Rate(views, rank = 2)
fit$status
fit$rate_multiplier
fit$rate_zero_certified
```

Remove `rate_c_g=1` from older launchers to inherit the new default, or set
`rate_c_g=0.05` explicitly. An old saved configuration containing 1 keeps 1.
The updated R and SLURM launchers in `inst/cluster/` require 0.2.14, remove that
old override, and retain the supplied signal of 0.8. They also use the existing
six-hour SGCA budget with infinite iteration caps. Use a fresh output directory
so completion markers from earlier runs cannot skip the new run. Copy the
updated launchers into the experiment directory expected by the SLURM scripts.

Further calibration is needed for materially different sample sizes, data
scales, numbers of views or population models. If sparse support recovery is
the primary target, the current results show that eliminating the zero estimate
is insufficient. A separately validated support-focused procedure would be a
different statistical choice, not an automatic fallback in this rate method.

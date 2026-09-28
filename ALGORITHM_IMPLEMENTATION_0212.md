# EGCAR 0.2.12: compressed ADMM and matrix-free spectral recovery

This release implements the proposed changes to Algorithms 1 and 2 of
`gca_scalability.pdf`, building on version 0.2.11. It preserves the two separate
statistical objectives, covariance divisor, ridge scaling, metric floor,
output thresholds, residual balancing, and cross-validation scoring rule.

## Changes

| Component | Implementation |
|---|---|
| L21 group and dual state | Store one endpoint array H=G+V and a row multiplier a, with G=aH and V=(1-a)H. |
| Algorithm 2, lines 5–9 | Form the C update from (2a-1)H, accumulate incident row norms during the C pass, then stream the endpoints once more for proximal updates and residuals. |
| Adaptive augmentation | For mu_new=f*mu_old, set b=a+(1-a)/f, H_new=b*H, a_new=a/b. This leaves G unchanged and rescales V by 1/f. |
| Legacy warm starts | Accept view-wide G/V or four endpoint arrays. Arbitrary legacy state undergoes one ordinary proximal update before compression; it is not projected to a different initial state. |
| Stopping rule | Preserve both endpoint primal constraints and the sum of endpoint auxiliary changes in the dual residual. Relative-tolerance norms and optional objective histories are preserved. |
| Support selection | Accumulate row norms from edge blocks; no full C is needed for support extraction. Existing output thresholding is preserved. |
| Leading eigenpairs | RSpectra receives a symmetric matrix-vector callback for S^(1/2) C S^(1/2), requesting exactly the r largest algebraic eigenpairs. |
| Metric square roots | Store covariance bases and scalar spectral weights. Apply both metric powers without materializing square-root matrices; retain the full covariance range and the ridge/floor on its complement. |
| Matrix-vector product | A registered C++ routine reads selected entries of each edge directly, without copying selected submatrices or assembling the full operator. |
| Failed partial solve | Retry with a larger Krylov space. If it still fails, report an error rather than compute a full EVD. A positivity threshold unresolved within the Ritz residual is reported as an invalid loading. |

The reference backend remains dense for independent diagnostics. Explicit
`partial_eigen=FALSE`, a user-selected larger `partial_eigen_min`, or a request
for all eigenpairs permits dense recovery. A 2x2 problem uses an analytic
leading eigenpair because RSpectra requires dimension at least three.

## Memory scope

Let q=sum(p_k*p_l) over k<l and p=sum(p_k). Main group ADMM state changes from
5q doubles (C plus two copies each of G and V) to 3q+p doubles (C plus two H
endpoint arrays and row multipliers). This is approximately 40% less state
memory. In the large benchmark, the retained R warm-start state was
96,802,984 bytes before and 58,104,864 bytes after, including R object overhead.

This is not a 40% reduction in total process memory. Prepared covariance
blocks, cross-covariances, spectral caches, reported coefficients, native/R
conversion copies, solver workspaces, and optional histories still occupy
memory. q remains quadratic for balanced views. Multiple independent CV
workers still have their own mutable states. Dense diagnostic output is
allocated if explicitly requested with `keep_full_C=TRUE`.

## Validation

The optional worker capability probe now calls the documented
`parallelly::supportsMulticore()` rather than an unavailable function in
`parallel`. `parallelly` is already a dependency of `future`; the existing
`EGCAR_MULTICORE=0` override is preserved. See
https://parallelly.futureverse.org/reference/supportsMulticore.html.

The source package was compiled and installed on Ubuntu 24.04, x86-64, with
R 4.3.3 and the system BLAS/LAPACK. The installed-package test suite passed
767 expectations with no failures or errors. One opt-in third-party comparison
suite was skipped; the separate serial/multisession and supported forked EGCAR CV tests ran.

`R CMD check --no-manual` completed with zero errors and zero warnings. Its
three informational notes concern unavailable optional suggested packages,
the explicit C++14 specification, and compiled-library size.

Tests include independent dense ADMM comparisons for tall and wide views,
zero columns, singleton dimensions, zero and positive penalties, adaptive mu,
arbitrary legacy starts, compressed starts, cross-backend continuation,
objective histories, unmodified inputs, selected-submatrix products,
metric-root actions, algebraic eigenvalue ordering, partial-solver failure,
requested ranks 1/2/5, and end-to-end experiment/CV behavior.

The tests also exposed and fixed a pre-existing failed-CV output bug: removing
a named NULL fit field let R partially match `$fit` to `fit_time`. Failed refits
now retain an explicit named NULL. Older fixtures were updated to request
retained CV rows explicitly and to compare identically preprocessed inputs.

## Measured timings

Each row reports the median of three runs using identical data, penalties,
fixed iteration caps, and serial compiled fitting. Preparation is excluded;
loading-factor construction is included. The smaller case uses n=50,
p=(320,360,400), r=2, 40 iterations. The larger case uses n=30,
p=(800,900,1000), r=2, 20 iterations. Both versions already use a partial
spectral solver on these dimensions; 0.2.11 first assembles its dense operator.

| Total p | Method | Total seconds, 0.2.11 | Total seconds, 0.2.12 | Speed ratio |
|---:|---|---:|---:|---:|
| 1080 | L11 | 1.430 | 1.400 | 1.02x |
| 1080 | L21 | 1.456 | 1.333 | 1.09x |
| 2700 | L11 | 4.100 | 2.197 | 1.87x |
| 2700 | L21 | 5.577 | 2.980 | 1.87x |

For p=2700, loading extraction decreased from 2.113 to 0.151 seconds for L11
and from 3.398 to 0.832 seconds for L21. ADMM solve time itself changed little
in these cases. The timing benefit depends on support size, covariance rank,
requested rank, spectral gaps, BLAS and hardware; these are not universal
speed guarantees or a benchmark of the user's full SLURM/CV sweep.

Across both benchmark sizes, selected sets agreed exactly. Maximum coefficient
disagreement was 1.1e-17, eigenvalue disagreement 1.6e-14, and eigenvector
projector disagreement 6.5e-14. Comparisons use projectors to allow eigenvector
signs and rotations within the recovered subspace.

Raw timing rows, numerical comparisons, test summaries, process peak-RSS
measurements and R session information are in `inst/validation/performance_0212_*`.
Peak RSS covers each entire benchmark R process, including preparation and
both penalty families; it is not per-fit native allocation accounting.

## Installation and compatibility

Install Rcpp, RcppArmadillo, RcppEigen and RSpectra, then install the supplied
`egcar_0.2.12.tar.gz` as an R source package. Restart R and CV workers afterwards.
The native state API marker is now 3, so stale compiled code is rejected.

Existing fitting calls continue to work. `egcar_control()` enables compressed
retained L21 state and matrix-free leading eigenpairs by default. Fits now
expose `solver$Hk`, `solver$Hl`, and `solver$a`. Pass
`egcar_control(compact_state=FALSE)` to obtain the old `solver$G`/`solver$V`
output; the optimized iterations still use compressed state. Existing old
fits and both old warm-start formats remain accepted. No rank-r restriction
is imposed on C and no approximate coefficient screening is introduced.

To repeat the benchmark with two installed versions:

```sh
Rscript tools/benchmark_0212.R /path/to/old_R_library /path/to/new_R_library results 3
Rscript tools/benchmark_0212.R /path/to/old_R_library /path/to/new_R_library results_large 3 large
```

To build and check the package:

```sh
Rscript tools/check.R /path/to/egcar
```

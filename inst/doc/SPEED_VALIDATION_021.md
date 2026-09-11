# EGCAR 0.2.1: speed-only validation report

Prepared 2026-09-10. Baseline: the supplied `egcar_0.2.0_github.zip`.

## Scope and source audit

All 27 existing files in `R/` are byte-for-byte identical to version 0.2.0.
`NAMESPACE`, all help pages, native registration, the native function signature,
and the existing standalone script are unchanged. No existing file was deleted.
The changes to existing implementation code are confined to the EGCAR native
solver in `src/egcar_native.cpp` and its standalone testing mirror. Build
metadata adds the header dependency RcppEigen. Installation/check helpers,
release documentation and additional tests/examples have been updated.

The L11-only and L21-only objectives, positivity of EGCAR CV grids, CV score,
fold sharing, warm-start conventions, stopping rules, numerical cutoffs,
loading extraction, SGCA/RGCCA/SGCCA/MultiCCA implementations, and oracle
calculations are retained. Floating-point reassociation can change final digits
or a numerical near-tie; mathematical equivalence is not bitwise equivalence.
The full simulation defaults were intentionally not changed in this revision.

## Native correctness tests executed

1. **120 independent dense-reference cases passed.** These include full and
   reduced covariance bases, a zero covariance view, unequal dimensions,
   zero/positive/large penalties, warm starts and adaptive augmentation.
   The maximum induced row-sum norm difference over all state matrices was
   8.8267584213985728e-11. The maximum objective discrepancy was
   3.5935698861067067e-12. The maximum reported residual difference was
   5.0098813986210189e-15. Iteration counts, convergence flags and final
   augmentation values matched in every case.
2. **864 comparisons with the exact version 0.2.0 native core passed.**
   Cases vary full/singular/zero-view/multiview inputs, warm starts,
   penalty family, coefficient, augmentation, adaptation and residual-check
   frequency. Maximum absolute elementwise state difference:
   1.0658141036401503e-14. Maximum residual/tolerance difference:
   3.6706748751669238e-15. Maximum history difference:
   4.801714581503802e-14. All stopping/adaptation decisions matched.
3. **8 native mapped-context cases passed.** Non-owning Armadillo views shared
   the source matrix storage, the solver produced matching results, and the
   source matrices were exactly unchanged. This tests native memory-view
   semantics, not Rcpp's R-object protection, R garbage collection or ALTREP.

The JSON reports and timing CSV are under `inst/validation/`. Reproducible native
source programs are in `tools/native/`. They use the same core text as the
package, not a separate reimplementation of its optimized arithmetic.

## Native-kernel timing experiment

Both versions used identical generated inputs, a penalty of 0.03, augmentation
0.7, 150 iterations, zero residual tolerances and disabled adaptation. This
fixes the computational workload rather than comparing differing stopping
iterations. Timings are medians over seven interleaved batches, alternating
old/new evaluation order. Each batch contains 30, 8 or 3 solves depending on
dimension. Reported times include native state copying and solver work, but
exclude covariance preparation, the Rcpp boundary, loading extraction, CV,
package loading and worker startup. These are NOT end-to-end R benchmarks.

For residual checks every five iterations, speed ratio = old time / new time:

| n | View dimensions | Total p | L11 ratio | L21 ratio |
|---:|:---|---:|---:|---:|
| 120 | 4, 4, 4 | 12 | 1.65x | 1.52x |
| 120 | 15, 15, 15 | 45 | 2.92x | 2.56x |
| 120 | 30, 30, 30 | 90 | 2.48x | 2.22x |
| 120 | 60, 60, 60 | 180 | 2.36x | 2.37x |
| 20 | 150, 150, 150 | 450 | 1.09x | 1.07x |

For example, at n=120 and three 15-variable views, the L11 kernel took
0.00535509 s previously versus 0.00183677 s after the update; the L21 kernel
took 0.00489881 s versus 0.00191062 s. Larger CV runtime reductions cannot be
inferred from these numbers because loading/scoring, R overhead and process
startup were not timed here. The p=450 scenario uses reduced bases (n=20).

Environment: Linux, g++ 14.2.0, C++14, -O3, Armadillo 14.2.3, Eigen 3.4.0,
Debian reference BLAS/LAPACK 3.12.1. BLAS/OpenMP thread environment variables
were set to one. Eigen parallelism is disabled in the core. Timing ratios are
hardware-, BLAS-, compiler-, dimension- and workload-dependent. They are not a
claim that the heuristic Eigen size cutoff is optimal on every machine.
R installation uses the user's RcppArmadillo/RcppEigen versions instead.

## Tests supplied but not executed

**R and Rscript are not installed in the available environment.** Accordingly,
none of the following has been executed here: package source installation,
compilation against Rcpp/RcppArmadillo/RcppEigen, native registration in R,
R-to-native context protection or input-immutability tests, R CMD check,
R-level EGCAR CV, comparison-package wrappers or multisession workers.
No claim is made that the full ten-method script has already succeeded.

`tests/testthat/test-speed-021.R` adds R checks for compiled-versus-R states,
input preservation, CV candidate scores and supplied folds. The existing test
suite is retained. The new installed example
`inst/examples/06_small_all_methods_check.R` attempts all ten methods on one
configuration (n=60, p_total=12, r=1, signal=0.8), with five folds and five
workers. It uses small tuning grids, not relaxed solver tolerances. It checks
configured equal worker allocation, records finite loading availability,
reports convergence separately and raises an error when a method lacks a
usable loading. The allocation file is configuration evidence, not CPU
utilization monitoring. This example should be run before a large study.

## Installation and distribution

Add RcppEigen to the existing source-build dependencies. A source build needs
a compatible R/C++ toolchain. No SMUT/EfficientCCA installation, third-party
headers, native executable or compiled package binary is included. The ZIP
contains one GitHub-ready `egcar/` source directory; the tar.gz is a source
archive assembled here, not the output of an executed R CMD build/check.
The existing author/maintainer placeholder must be replaced for publication.

See `EFFICIENTCCA_REVIEW.md` for the inspected upstream sources, adopted
computational ideas and deliberately unimported statistical/solver changes.

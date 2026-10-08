# EGCAR 0.2.17 validation and limitations

Validated on 2026-10-08 using local R 4.3.3, a single BLAS thread, and the released 0.2.16 sources as the baseline. The cluster workflow defaults to R/4.2.0 and runs an install/preflight check there before starting experiments. No Slurm jobs were submitted from this environment.

## Model and scope

Only the L21 fitting path changes. The new objective is the directed, zero-diagonal-block quadratic loss with penalty `lambda * sum_i ||C[i,]||_2`. The sparse ADMM copy is averaged with its transpose once after the solve, before loading extraction. The loss has half weight per directed pair so that it agrees with the former loss on symmetric C.

The old sum of view-group norms already equals the global row norm when C is symmetric. Delaying symmetry therefore changes the feasible set and the estimator. Convergence describes the directed problem, not optimality of the final symmetric projection for the old constrained problem. Symmetrization can expand the set of nonzero rows. The retained rate multiplier 0.05 has not been recalibrated for this estimator.

L11, SGCA and the other benchmark algorithms, the data generator, the seed rules, error metrics and loading extraction were not changed. No other method is executed by the accompanying cluster workflow. Original L21 warm-start states are rejected to avoid mixing formulations.

## Numerical and workflow checks that passed

- 51 selected package tests, 978 assertions, zero failures, errors, skips or warnings; the exact test list and counts are in package_test_results.csv.
- Native solver agreement with an independent explicit R implementation on unequal blocks, full-rank and rank-deficient/wide designs, and five-block designs; checked adaptive penalty updates, iteration histories and warm starts.
- Independent identity-Hessian closed form: row thresholding at twice lambda, followed by the final symmetry projection. This checks the loss and penalty factors independently of the ADMM recurrence.
- Exact returned symmetry, zero diagonal blocks, compact state and correct full-space nullspace contributions.
- Existing L11 reference regressions, ranks 1/2/5 loading tests, and timing tests that include loading extraction. SGCA paper-specific in-process regressions also passed.
- All 2,850 dataset mappings, including population, observation and fold seeds and both original method IDs 5/6, agree with the original seed rules. New tiers extend those rules to signals 0.3 and 0.5.
- An original tier dataset and folds were reproduced exactly; immutable original-bundle references, regenerated-data agreement with original metric hashes, idempotent preparation and tamper rejection passed.
- Local rate and five-worker/five-fold CV preflight fits both returned status ok with valid loading vectors and consistent timing. All six experiment grids loaded.
- Mock Slurm checks passed for six groups, separate rate/CV chains, bounded arrays, CPU allocation and aggregation dependencies. No real scheduler execution was simulated as a completed cluster experiment.
- Aggregation checks passed for committed result markers, package version, all expected cells (including cells with no results), usable-result counts and Sigma0-error/runtime plots. Shell syntax and Python compilation passed.

The full historical suite was also attempted, but unrelated SGCA subprocess-supervision tests encountered local callr/ps/processx process-information failures. They are not included in the passing-test claim above. Optional RGCCA/PMA numerical paths were not retested. Full-scale 5,700-fit experiments have not been run here.

## Runtime and memory measurements

The user's requirement of at least the previous speed and peak-memory efficiency is **not established; the tested new solver was slower**. These measurements are disclosed so this package can be evaluated as a changed statistical estimator rather than represented as a speed upgrade.

Three centered Gaussian views, n=150, rank=1, seed=732; `lambda=0.05*sqrt((2*p+log(3*p))/n)`, max_iter=2000 and check_every=5. Each entry is the median of three solver-plus-loading times in one process; preparation and data generation are excluded from that time. Each repetition returned valid loadings and the reported fits converged. The legacy and new implementations used the same prepared data and installed libraries. This is a small synthetic benchmark, not the Toeplitz production experiment.

| Variables per view | Legacy time (s) | New time (s) | Legacy iterations | New iterations | Legacy process peak RSS (KiB) | New process peak RSS (KiB) |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.633 | 1.213 | 50 | 70 | 223300 | 224272 |
| 500 | 4.783 | 14.431 | 15 | 30 | 263168 | 276192 |

Peak RSS includes the whole R process, package, prepared covariance caches and all three repetitions; it is not an allocation count. It was slightly higher for the new solver in these runs. CPU timing and RSS vary by hardware and runtime. Two independent directed updates and different iteration counts help explain the slowdown. Since the estimators differ, these numbers do not establish equal statistical accuracy.

The native persistent ADMM state decreases from `3*q+P` to `2*q+P` doubles, with `q=sum_{k<l} p_k*p_l` and `P=sum_k p_k`. This is a state-storage reduction, not a claim about total peak RSS. Covariance caches, scratch, input objects and returned results are additional. Final symmetry uses block operations without constructing a dense P-by-P matrix. Thin spectral bases, cached denominators, partial loading EVD, warm starts, single-threaded BLAS inside workers, and parallel CV folds are retained. An exact zero-solution KKT check skips iterative solves when lambda already guarantees zero; the resulting unavailable loadings remain recorded as failures.

For a portable timing rerun after installation, use `Rscript validation_report/benchmark_l21.R old 100 3 150` and the same command with `new`, optionally wrapped in `/usr/bin/time -v` on Linux. The archived measurements used a resource.getrusage wrapper for peak RSS. Broader convergence-matched performance testing would be required for a general speed claim.

## Experiment coverage and comparisons

The supplied original tier scripts covered ranks 1–5 and ten repetitions at signal 0.8 only. The supplied original consistency archive covered all three signals 0.3/0.5/0.8 and ranks 1–5, with ten repetitions. These are code-level coverage findings, not statements that remote jobs completed successfully.

The new workflow has 600 tier datasets and 2,250 consistency datasets, each fitted separately by EGCAR-L21-rate and EGCAR-L21-CV: 2,850 datasets and 5,700 fits in total. CV parallelizes over five folds; each fold evaluates its parameter candidates sequentially. Fit times include the final symmetry step and loading-vector calculation; CV total time also includes tuning.

At submission, original source paths and run IDs are frozen. Existing input bundles are checked and referenced without copying their data. If absent, the original seeds regenerate them; available old metric hashes must match. If no old hashes exist, identity is explicitly recorded as seed-based and unverified against an old byte-level dataset. Same seeds alone are not a guarantee across changed R/RNG/numerical environments.

Signals 0.3/0.5 are new in the tiers: no old results for the other methods exist at those settings in the supplied tier design. This L21-only rerun does not supply those missing comparison fits. Old and new L21 outputs are stored separately. Successful-fit plots exclude numerical failures and report the usable denominator; completion is not treated as proof of numerical validity.

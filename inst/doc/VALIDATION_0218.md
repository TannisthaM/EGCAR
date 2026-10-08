# EGCAR 0.2.18: validation and performance limits

Validated locally on 2026-10-08 with R 4.3.3. The cluster bundle uses R/4.2.0 by default and runs installation and a preflight there. No real Slurm jobs were submitted from this environment.

## What changed

The statistical estimator is the directed global-row-L21 model introduced in 0.2.17. Its penalty is lambda times the sum of row Euclidean norms of C; zero diagonal blocks are enforced. After ADMM, the sparse proximal estimate Z is projected to (Z+Z')/2 and loading vectors are computed. The former view-group sum already equalled global row L21 when C was symmetric. Delaying symmetry changes the feasible set and hence the estimator; final symmetry does not turn this into the old symmetry-constrained optimization problem. Projection can fill previously zero rows.

0.2.18 optimizes that directed problem using right spectral coordinates and an exact compressed representation of complementary components. It uses ADMM over-relaxation alpha=1.6 and relative-residual balancing. The original mu argument remains the initial numerical penalty. The statistical coefficient lambda, full-space absolute/relative stopping tolerances and iteration caps are retained. The stopping check remains on unrelaxed primal feasibility and the dual residual, not on a smaller surrogate dimension.

Only L21 fitting and its warm-state handling change. The data generator, folds, L11 numerical solver, SGCA/RGCCA/SGCCA/MultiCCA algorithms, error metrics and loading extraction remain unchanged. The previous L21 rate multiplier 0.05 is retained for comparability; it has NOT been recalibrated for the new estimator. Every recorded method result identifies package 0.2.18, the formulation, alpha and the mu-balancing rule.

## Checks passed

- 53 selected package tests, 1,016 assertions, zero failures/errors/skips/warnings. See package_test_results.csv for the exact list. This is not a claim that the entire historical suite or R CMD check passed.
- Native agreement with an independently implemented full-matrix R ADMM solver: coefficients, full-space histories/residuals, adaptive mu, unequal blocks, rank-deficient/wide cases, five-block cases and cross-backend warm starts.
- Explicit nonzero covariance-nullspace forcing, including lambda=0, fixed and adaptive mu, and warm starts; the complementary terms are retained rather than dropped.
- An independent identity-Hessian closed form (row threshold at 2*lambda, then symmetry), zero-solution certificates and exact returned symmetry.
- L11 reference regressions, requested loading ranks 1/2/5, partial-versus-dense loading eigensolves and full-loading timing checks.
- SGCA paper-specific in-process regressions passed. Optional RGCCA/PMA numerical paths and SGCA subprocess-supervision tests were not rerun; the latter had local callr/process-information failures in the prior validation attempt. Their source code was not changed.
- All 2,850 dataset mappings and both method seeds checked against the original seed rules. The dataset generator and benchmark source files were also compared byte-for-byte with 0.2.16.
- Original tier data/folds reproduced, original bundles referenced without copying, original metric hashes checked when available, repeated preparation was idempotent and modified source bundles were rejected.
- Local rate and five-worker/five-fold CV preflight returned status ok, valid loadings and consistent timing; all six grids loaded. Five distinct worker process IDs were observed. A separate ten-block/rank-five fit also converged with all five valid loading vectors.
- Mock scheduler checks passed for bounded arrays, independent L21-only queues, five-worker CPU allocation and aggregation dependencies.
- Aggregation checks passed for commit markers, versions, input identities, all expected cells including empty cells, usable counts, and error/runtime plots. Combined tier plots include p=10,100,500,1000. Shell/Python syntax checks passed.

## Runtime and whole-process memory measurements

The optimized solver improves the earlier 0.2.17 implementation, but a strict guarantee of never exceeding 0.2.16 runtime or total peak RSS is NOT established. Some measured settings still have slightly higher whole-process memory, and the cold-start p=500 Toeplitz check was roughly tied in runtime. Do not interpret the state-storage reduction as a universal peak-memory or speed guarantee.

The first table uses three Gaussian views, n=150, rank=1, seed=732, lambda=0.05*sqrt((2*p+log(3*p))/n), max_iter=2000 and check_every=5. Each runtime is the median of three solver-plus-loading fits in one process. The actual 0.2.16 and 0.2.18 installations were used. Both versions had RSpectra loaded before this measurement; data generation and preparation are excluded from time. Every fitted result converged and returned valid loadings. Peak RSS includes the whole process and all three repetitions, including setup. These are small diagnostics, not completed production experiments.

| p per view | 0.2.16 time (s) | 0.2.18 time (s) | Old/new iterations | 0.2.16 peak RSS (KiB) | 0.2.18 peak RSS (KiB) |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.536 | 0.386 | 50 / 35 | 203920 | 204456 |
| 500 | 4.518 | 3.803 | 15 / 20 | 292512 | 265752 |
| 1000 | 11.504 | 8.906 | 10 / 15 | 482828 | 388920 |

This was about 28%, 16% and 23% less time, respectively. Peak RSS was 0.3% higher at p=100 and 9%/19% lower at p=500/1000. Small process-memory differences can depend on R allocation and garbage collection; no confidence interval is inferred from this table.

The next table uses the ACTUAL original tier seed rule, rank=1, replication=1, n=150, signal=0.8, Toeplitz parameters (0.5,0.7,0.9), and five active variables per block. Each entry is a single fresh-process fit and INCLUDES first-use loading-library initialization inside the fit timer. This preserves the original method-job timing boundary. Dataset SHA256 hashes agree exactly between versions and are included in the logs. Both versions converged with valid loadings. Population generation and preparation remain outside the fit timer but inside process peak RSS.

| p per view | 0.2.16 full fit (s) | 0.2.18 full fit (s) | Old/new iterations | 0.2.16 peak RSS (KiB) | 0.2.18 peak RSS (KiB) |
|---:|---:|---:|---:|---:|---:|
| 100 | 1.158 | 1.022 | 70 / 40 | 181068 | 185328 |
| 500 | 6.345 | 6.409 | 20 / 40 | 296960 | 308152 |

The single p=500 cold run was about 1% slower and used about 4% more peak RSS. The p=100 run was faster but used about 2% more RSS. These measurements are retained rather than selectively reporting only improvements. Full CV performance, all ranks/signals and all sample sizes have not been benchmarked against 0.2.16. The changed estimators need not have equal statistical accuracy just because both numerical stop tests pass.

## Storage and arithmetic

Let q=sum_{k<l}p_k*p_l, P=sum_k p_k, and d_k=ncol(Q_k). The old native persistent ADMM iterates used 3q+P doubles. New iterates use s+2P, where s=sum_{k<l}(p_k*d_l+p_l*d_k). The extra P-vector stores the complete right-nullspace contribution, not an approximation. For equal p and d, the dense-state reduction is substantial when d is much smaller than p.

Covariance caches, right-hand-side remainder caches, work buffers, R objects and the returned symmetric blocks are additional. Near-full-rank cases can have sizable remainder caches; the iterate formula alone must not be used as a total-memory bound. Cold arrays are allocated directly in C++. When a final fit does not retain warm state, each directed edge is released as its symmetric result is produced. No dense P-by-P C is assembled.

For equal-sized blocks, the leading dense multiplies per unordered edge change from approximately 2*(p^2*d+p*d^2) to 4*p*d^2 per ADMM iteration. Iteration counts can change, so this is not a total-runtime guarantee. There is additional one-time projection/reconstruction work per fit. Loading extraction still uses the existing matrix-free partial EVD and is included in all fit times.

Implementation notes and derivation are in L21_IMPLEMENTATION.md. General background: Boyd et al., ADMM (https://stanford.edu/~boyd/papers/admm_distr_stats.html), and Wohlberg, residual balancing (https://arxiv.org/abs/1704.06209). The specific compressed representation and implementation are derived for this package and tested against the full-space solver.

## Coverage and seeds

The original tiers included ranks 1–5 and ten repetitions at signal 0.8 only. Original consistency included ranks 1–5, all signals 0.3/0.5/0.8 and ten repetitions for each requested setting. This verifies scripts, not completion of remote jobs.

The new workflow has 600 tier datasets and 2,250 consistency datasets, with independent rate/CV jobs for each: 2,850 datasets and 5,700 fits. Seeds and original method IDs 5/6 are retained. Existing original bundles are verified and referenced directly. Missing inputs are regenerated from the original seed rules, with old metric hashes checked when available; otherwise the result explicitly says that old byte-level identity was unavailable. Same seeds alone do not certify equality across changed RNG/numerical environments.

The new 0.3/0.5 tier settings have no old other-method results in the supplied tier design. Their new L21-only fits do not create those missing benchmark comparisons. The 5,700 production fits have not been run locally. Failed or nonconverged fits remain explicit failures; successful-fit averages display their usable denominator.

# EGCAR 0.2.18: global row-L21, followed by symmetrization

This folder runs only EGCAR-L21-rate and EGCAR-L21-CV (original method IDs
5 and 6). Each method has an independent Slurm queue. Existing experiments
and their R libraries are not modified.

## What the earlier archives actually covered

| Study | Variables per block | Sample sizes | Ranks | Signals | Replications |
|---|---|---|---|---|---|
| Original tier 1 | 10, 100 | 150 | 1–5 | **0.8 only** | 10 |
| Original tier 2 | 500 | 150 | 1–5 | **0.8 only** | 10 |
| Original tier 3 | 1000 | 150 | 1–5 | **0.8 only** | 10 |
| Original consistency K=3,5 | 100 | 10,100,500,1000,10000 | 1–5 | 0.3,0.5,0.8 | 10 |
| Original consistency K=10 | 100 | 100,500,1000,10000,50000 | 1–5 | 0.3,0.5,0.8 | 10 |

These statements verify the supplied scripts, not the completion of remote jobs.
The new tiers add signals 0.3 and 0.5. There are no corresponding old tier
benchmark outputs at those signals. This L21-only workflow cannot create the
missing results for the other eight methods.

## New estimator

On a directed coefficient matrix C with zero diagonal blocks, minimize

    (1/4) sum_{k != l} tr(C_kl' S_kk C_kl S_ll)
      - (1/2) sum_{k != l} <S_kl, C_kl>
      + lambda sum_i ||C[i, ]||_2.

The half weighting per directed pair preserves the old loss scale when C is
symmetric. Lambda is still the CV-selected coefficient multiplying the norm.
The row penalty is used; row and column penalties need not be equal before
symmetrization. ADMM fits the directed problem. Its sparse proximal
copy Z is returned as (Z + Z') / 2, stored as one block per unordered pair.
Loadings are then extracted using the existing partial eigensolver.

In the old estimator, sum_k ||M_k(C)||_21 already equalled ||C||_21 on the
symmetric zero-diagonal-block matrix. The substantive change is relaxing
symmetry during fitting and imposing it afterward. This is a different
estimator. The final projection can fill rows that were zero in Z, and its
output is not claimed to minimize the old symmetry-constrained objective.
ADMM convergence is a certificate for the directed problem before projection.

Only L21 fitting changes. Simulation, folds, L11, SGCA, RGCCA, SGCCA, MultiCCA,
the error metrics, and loading extraction remain unchanged. The existing
L21 rate multiplier 0.05 is retained for comparability; it was not recalibrated
for this new estimator. The existing CV grid, tolerances, and iteration caps
are retained. Numerical failures remain recorded failures.

## Memory and runtime

The native L21 solver works in right spectral coordinates. An outgoing block
has p_k by d_l entries, with d_l the retained covariance-basis dimension.
Orthogonal right rotations preserve row norms. The entire complementary
component is also retained: for global row i in block k it is a scalar h_ki
multiplying the corresponding row of S_kl(I-Q_l Q_l'). It is never silently
screened out. One scalar per global row suffices because the same row proximal
multiplier acts on every outgoing block.

Persistent ADMM iterate storage is s+2P doubles, with
s=sum_{k<l}(p_k*d_l+p_l*d_k), P=sum_k p_k. The old native solver kept
3q+P, q=sum_{k<l}p_k*p_l. Covariance caches, right-hand-side remainders, scratch,
R objects and returned C are additional; this state formula is not peak RSS.
Full C is never assembled during ADMM or for final symmetrization.
Cold-start arrays are allocated natively. For fits without retained state,
edges are released as their final symmetric blocks are reconstructed.

The new solver also uses over-relaxation and relative-residual balancing of
ADMM's numerical mu. The coefficient lambda, the statistical loss, full-space
stopping tolerances and iteration caps are unchanged from the directed model
in 0.2.17. These numerical accelerations change the iteration trajectory.

Read VALIDATION.md for actual comparisons with 0.2.16, including limitations.
Neither fewer state arrays nor faster measured cases proves that every future
setting will have lower runtime and lower whole-process peak memory.

Thin spectral bases, partial loading EVD and CV warm starts are retained.
CV uses five parallel fold workers, with candidates sequential within each
fold. The fit timer includes ADMM, final symmetrization and loading vectors;
CV total_time also includes tuning. The method jobs retain the original
first-fit dependency-loading behavior. No setup cost is moved out of the timer
merely to improve a comparison.

## Inputs, seeds, and pairing

Original configs are preserved verbatim in r/original_*_config.R. The new
workflow retains original dataset IDs separately from its own scheduling IDs.
It never uses Slurm IDs or job order as random seeds.

At submission it records the original data/result paths in
runs/<run>/source_roots.csv. Available original bundles are referenced directly,
without copying or modifying them; keep those source bundles in place.
Checksums, dimensions, signals, data hashes, fold hashes, and seeds are checked.
If an input bundle is unavailable, the dataset is regenerated from the original
seed rule. Where an old metric is available, its hashes must match exactly.
Otherwise the result is honestly marked same_seeds_original_hash_unavailable.

For tiers, the original rule was

    task = 20260907 + 100000*rep + 1000*r
    population = task + 100000 + 1000*r
    data = population + 77
    folds = task + 1000000 + 10000*r + n
    method = folds + 1000*original_method_id.

The new 0.3/0.5 tier settings use these same seed rules: common random numbers
across signals, ten independent replications within each setting. A seed alone
does not mean that changing signal produces the same observations.

Consistency uses its exact original ordering: signal, n, rank, replication.
Each K has 750 datasets. With global_id=(K_index-1)*750+old_dataset_id,

    task = 20261002 + 1000000 + 100*global_id
    population = 20261002 + 10000*K_index + 100*r + rep
    data = task + 1; folds = task + 2
    method = task + 10 + original_method_id.

If the original locations differ, set EGCAR_OLD_TIERS_ROOT,
EGCAR_OLD_TIERS_RUN, EGCAR_OLD_CONSISTENCY_ROOT, and/or
EGCAR_OLD_CONSISTENCY_RUN before submission. The tier run defaults to
separate_20260929T064646. Consistency reads its old LAST_RUN.txt once at
submission. Audit the printed original run names and source_roots.csv.

## Upload and launch

On the Mac, after downloading egcar_l21_0218.zip:

```bash
cd "$HOME/Downloads"
scp egcar_l21_0218.zip tannistha@midway3.rcc.uchicago.edu:
ssh tannistha@midway3.rcc.uchicago.edu
```

On the cluster:

```bash
unzip -n "$HOME/egcar_l21_0218.zip" -d "$HOME"
cd "$HOME/egcar_l21_0218"
sha256sum -c CHECKSUMS.sha256
bash plan.sh
bash download_dependencies.sh
bash submit_all.sh
bash status.sh
```

All 6 groups are launched. Initial submissions are an install/preflight job
plus six preparation launchers. Each group prepares bounded batches, then
launches independent method-5 and method-6 chains. Future batches are submitted
as earlier batches end; 5,700 jobs are not all placed in the queue at once.
Job IDs are retained automatically in runs/<run>/jobs.tsv. Prefix: eg218l_.

| Group ID | Study | Datasets | L21 tasks |
|---|---|---:|---:|
| 1 | tier 1, all three signals | 300 | 600 |
| 2 | tier 2, all three signals | 150 | 300 |
| 3 | tier 3, all three signals | 150 | 300 |
| 4 | consistency K=3 | 750 | 1500 |
| 5 | consistency K=5 | 750 | 1500 |
| 6 | consistency K=10 | 750 | 1500 |
| Total | | 2850 | 5700 |

Defaults: account pi-cdonnat, partition caslake, R/4.2.0; six CPUs for CV and
one for rate fits; one BLAS thread. Memory per fitting task is 20/40/50/16/24/64
GB for groups 1–6. Both method chains default to two concurrent tasks per group;
preparation defaults to three. Set EGCAR_METHOD_CONCURRENCY=1 and
EGCAR_PREP_CONCURRENCY=1 before submitting to reduce concurrent resources.
There is no automatic change to estimator rank, ridge, tolerance or penalty.

## Monitor, aggregate, download

```bash
bash status.sh --accounting
bash status.sh --detail
bash aggregate_now.sh
bash pack_results.sh
```

aggregate_now.sh queues aggregation jobs; wait for those jobs to finish before
packing their outputs. Aggregation is also queued automatically when both method queues for a group
end. Inspect completion_counts.csv: recorded outcomes include unsuccessful
fits. Plotted averages use only valid, converged estimates and show usable n/10.
Outputs include all_metrics.csv, failure records, input identities, CV diagnostics,
and Sigma0-error/runtime PDFs by signal and rank. Combined tier curves across
p=10,100,500,1000 are in results/<run>/tiers_combined/plots; consistency groups
produce curves against n. Group folders group1–group6
refer to group IDs in the table, not the actual number of data blocks.
The download archive excludes original or regenerated input bundles.

On the Mac:

```bash
cd "$HOME/Downloads"
scp tannistha@midway3.rcc.uchicago.edu:egcar_l21_0218/egcar_l21_results_latest.tar.gz .
mkdir -p egcar_l21_results
tar -xzf egcar_l21_results_latest.tar.gz -C egcar_l21_results
open egcar_l21_results
```

Use resume_submission.sh for interrupted submissions. retry_missing.sh GROUP
METHOD reruns only missing outputs (METHOD=0 prepares missing input, 5 is rate,
6 is CV); it never overwrites a committed numerical failure or a running chain.
Original results and methods are not rerun by these scripts.

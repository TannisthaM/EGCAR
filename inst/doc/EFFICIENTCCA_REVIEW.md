# EfficientCCA review and the EGCAR 0.2.1 speed-only update

## Sources inspected

Inspected the public `main` branch on 2026-09-10:

- https://github.com/ZixuanWu1/EfficientCCA
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/EfficientCCA.R
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/EfficientCCA_update.R
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/GCA/gca_to_cca.R
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/GCA/init_process.R
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/GCA/adaptive_lasso.R
- https://raw.githubusercontent.com/ZixuanWu1/EfficientCCA/main/R/GCA/utils.R
- https://raw.githubusercontent.com/cran/SMUT/master/src/matrixproduct.cpp

The linked project is a repository of research scripts, not an additional
package dependency in this release. No upstream implementation was copied.
The methods and benchmark algorithms remain the EGCAR 0.2.0 implementations.

## Findings from the sources

`EfficientCCA_update.R` caches covariance eigendecompositions, a spectral
cross-covariance and augmented denominators before traversing penalty values.
It calls the compiled `SMUT::eigenMapMatMult` operation, requests only the
leading singular vectors using `svds`, and parallelizes validation folds.
SMUT's multiplication implementation maps its input matrices into Eigen.

EGCAR 0.2.0 already cached the corresponding spectral information, used
warm-started coefficient paths, supported checked partial eigendecompositions,
and shared a fold-level worker pool. It additionally handled wide views with
thin covariance bases and the complete null-space correction. Replacing that
implementation wholesale was therefore not justified.

## Additional EGCAR-specific implementation changes

1. Mapped Eigen products now handle small-to-moderate dense products in the
   compiled ADMM projection and reconstruction. The maximum of the three
   multiplication dimensions must be greater than 8 and at most 64. Smaller
   products use Armadillo's tiny kernels; larger products retain the existing
   Armadillo/BLAS multiplication. These are dispatch heuristics, not a claim
   that one backend is fastest on every platform. Eigen's internal parallelism
   is explicitly disabled, so it does not create another worker pool.
2. Immutable double-valued covariance/spectral context matrices are mapped
   into Armadillo for the lifetime of the native call rather than copied on
   each CV candidate. Warm-start states remain owning copies and are never
   mutated in R. No pointer is returned, cached across calls, or exported to
   a parallel worker. Non-double context inputs retain the owning conversion.
3. A single reusable projection/reconstruction workspace is shared by the
   sequential edge updates. Its buffers grow as necessary and are reused
   throughout the solve. Denominators still change only when augmentation
   changes. Group-norm buffers are also reused.
4. Entrywise/group proximal operations, dual updates and full-space residual
   reductions are fused into fewer matrix traversals. The original residual
   checks, adaptation schedule and stopping tolerances are retained.

The small-matrix Eigen branch and mapping avoid an SMUT dependency. The only
new build dependency is `RcppEigen`, which supplies Eigen headers. The
existing `RcppArmadillo` dependency remains. Both separate penalty families
benefit for fixed, rate-scaled and CV fitting with `backend = "cpp"`.
The optimized R and dense reference backends are unchanged.

## Deliberately not imported

The upstream coefficient-change stopping rule is not substituted for EGCAR's
primal/dual residual rules. The two-view validation MSE is not substituted
for the multiview held-out score. Upstream group weights are not substituted
for the concatenated incident-row penalty, and no automatic cluster is created
at package load. The EGCAR L21 edge shift remains `2 * mu`. Both EGCAR CV grids
remain strictly positive. No rank constraint, screening, altered covariance
regularization, or new statistical estimator has been introduced.

All `R/` implementation files, `NAMESPACE`, `man/` help pages, original solver
controls, benchmark CV procedures, simulation routines and oracle definitions
are byte-for-byte identical to version 0.2.0. Only the EGCAR native source,
build-dependency metadata, installation/check helpers, release documentation
and added validation/example files change. The package's full-study defaults
are intentionally not changed as part of this speed-only revision.

## Small end-to-end check

Run `inst/examples/06_small_all_methods_check.R` after installing the package.
It requests one replication at n=60, p_total=12 (three views of size four),
r=1 and signal=0.8, with two candidate values per tuning dimension and five
shared folds/workers. It attempts all ten existing methods and records
finite-loading availability separately from solver convergence. Iteration
limits and tolerances are not shortened. Plotting is disabled for this check.
See `SPEED_VALIDATION_021.md` for checks executed during preparation and
limitations of the available runtime.

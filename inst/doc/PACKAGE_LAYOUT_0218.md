# Repository and documentation refresh, 2026-10-08

Based on the public `TannisthaM/EGCAR` main-branch commit
`c64e513d74d0f090cb18092234d5faca54f3d88a`, which already contained version 0.2.18.
The package version is retained because this update does not change the
numerical algorithms, public function signatures, defaults, or runtime accounting.

## Changes

- Replaced the accumulated release-history README with current installation,
  data requirements, method descriptions, small examples, and timing guidance.
- Added installed `quickstart.R` and `compare_methods.R` examples using public APIs.
- Added an RStudio package project, source/build ignore rules, GPL-3 license text,
  and contributor instructions. Removed a tracked macOS metadata file.
- Moved the historical algorithm and optimization reports into `inst/doc/`.
- Added repository/issue URLs to `DESCRIPTION`. Replaced the placeholder
  maintainer with the name and email in the repository's latest public commit;
  retained existing author/copyright attribution and the bundled MIT notice.
- Removed generated `Author`, `Maintainer`, and `Packaged` fields from source
  metadata; R derives author/maintainer information from `Authors@R` on build.
- Corrected stale package-level documentation and the placeholder GitHub installer.

## Checks performed

On Linux with R 4.3.3:

- Source build and installation succeeded, including native compilation.
- `R CMD check --no-manual --no-tests` with optional Suggests nonmandatory:
  **0 errors, 0 warnings, 4 notes**. Notes concern unavailable optional packages,
  the C++14 specification, installed library size, and a retained historical
  LaTeX source in `inst/doc/`. Help-page examples passed.
- All six fits in the installed quick-start script returned `status="ok"`
  and `converged=TRUE`, with valid loading extraction and prediction.
- The README's four-fit example and three-worker CV example ran successfully.
- The comparison example's two oracle calls, MultiCCA call (PMA 1.2-4), and
  finite-budget SGCA call ran successfully on its supplied small dataset.
  The supervisor emitted process-stat warnings in this environment but returned
  a valid SGCA result. Completion is not a claim of numerical convergence.
- README/example-index local links were checked.

The full historical regression suite was not rerun for this documentation-only
update. RGCCA and SGCCA were not executed: the optional RGCCA dependency chain
was unavailable in this R environment. No new performance benchmark or cluster
experiment was run. Existing algorithm validation is in `VALIDATION_0218.md`.

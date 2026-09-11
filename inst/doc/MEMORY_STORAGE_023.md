# EGCAR 0.2.3 memory/storage revision

The experiment defaults now prioritize low retained memory and compact disk output while preserving the statistical computations.

Defaults: `save_fits = FALSE`, `save_cv_fold_results = FALSE`, `save_loading_data = FALSE`, `save_compact_loadings = TRUE`, `retain_benchmark_fits = FALSE`, `make_loading_plots = FALSE`, and `loading_factor_cache_max = 4L`. Per-configuration loading PDFs can be regenerated later from the compact loading archives.

`compact_loadings/` stores only the raw loading matrices, truth, view sizes, active-row indices, identifiers, and method status. Full solver/package fits and populations are omitted from the main checkpoint unless `save_fits=TRUE`. Fold-level diagnostics are omitted unless `save_cv_fold_results=TRUE`. Full loading-plot RDS/CSV files are omitted unless `save_loading_data=TRUE`.

ADMM warm-start state is still kept within a CV path because the next penalty candidate needs it. Rate fits and selected full-sample refits use compact results after `C_hat` is formed. The public `egcar_fit()` API keeps its prior full-state behavior.

The SGCA initializer cache keeps the derived rank-r initializer and scalar convergence diagnostics rather than full Pi/H/Gamma matrices. Full benchmark fit objects remain available to public comparison calls, while the experiment runner defaults to not retaining them.

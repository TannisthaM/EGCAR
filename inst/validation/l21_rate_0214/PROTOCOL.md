# Prespecified L21 rate calibration for egcar 0.2.14

Written before running the calibration candidates. The supplied all_metrics CSV
and its ground-truth errors are not inputs to this study.

- Existing objective, solver, covariance divisor, row thresholds and rank are unchanged.
- Rate: c * sqrt((max_k(p - p_k) + log(p)) / n).
- Candidates: 0.05, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75, 1.00.
- Calibration: 36 independent simulations, crossing p/block = 10, 100, 500,
  1000; ranks = 1, 3, 5; signals = 0.3, 0.5, 0.8. Each has 150 training
  observations and 600 independent validation observations from its population.
- Three views, five active variables per view, Toeplitz correlations 0.5, 0.7,
  0.9, matching the experiment model. Seeds start at 100214000 and are explicitly
  listed by l21_pilot_design(). No original experiment seeds are reused.
- Select among candidates producing a valid, converged requested-rank loading
  in EVERY calibration scenario. Maximize the mean held-out score divided by
  rank, with equal weight per scenario. Exact ties favor the larger multiplier.
  No population subspace error or true support is used for selection.
- A strictly greater penalty than half the largest incident cross-covariance
  row norm certifies a unique zero solution. Such candidates are recorded as
  invalid_loading without running ADMM; they are not dropped from counts.
- Freeze the selection in selected_multiplier.txt before confirmation.
- Confirmation: 60 independent new populations, crossing the same dimensions
  and signals with ranks 1:5. Seeds start at 200214000. Evaluate only the frozen
  multiplier and the old multiplier 1. Report validity, support size, held-out
  score and population-metric subspace errors without revising the selection.
- This is a finite pilot for the supplied simulation regimes, not an optimality
  theorem or a guarantee for arbitrary data, dimensions or sample sizes.

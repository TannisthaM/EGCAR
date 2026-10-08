# Benchmark reproduction

The logs are the actual recorded local R 4.3.3 measurements. On Midway, load the same R module for both versions and point the portable scripts at each installed library. Output reports the package version actually loaded. Run sequentially on an allocated compute node with one BLAS thread; CPU/memory comparisons should use equivalent resources.

Examples, from the bundle root after `source env.sh`:

```bash
python3 validation_report/measure.py "$HOME/egcar_separate_0216/Rlib" 500 3 150
python3 validation_report/measure.py "$PWD/Rlib" 500 3 150
python3 validation_report/measure_tier.py "$HOME/egcar_separate_0216/Rlib" 500
python3 validation_report/measure_tier.py "$PWD/Rlib" 500
```

`bench.R` preloads RSpectra and reports the median of three Gaussian fits. `benchmark_tier.R` runs one original-seed Toeplitz fit, including first-use loading-library initialization. Both time through loading extraction. Neither includes data generation or preparation in fit time; process peak RSS includes them. Reading VALIDATION.md is necessary to interpret the different timing boundaries. These comparisons do not alter production scripts or submit new jobs.

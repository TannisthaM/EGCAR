# Working on EGCAR

The repository root is the R package root. Open `egcar.Rproj` in RStudio,
or work directly with R and a C++14 toolchain.

## Install and try the package

```r
install.packages("remotes")
remotes::install_local(".", dependencies = NA, upgrade = "never")
source(system.file("examples", "quickstart.R", package = "egcar"))
```

Use a fresh R session after rebuilding the native library.

## Check a change

```r
install.packages(c("testthat", "future", "future.apply"))
```

From the repository root:

```bash
Rscript tools/check.R .
```

This parses source/help files, builds the source package, and runs
`R CMD check --no-manual`. Read the resulting check log for warnings and
notes as well as errors. To include the optional external-method checks,
install `RGCCA` and `PMA` and set `EGCAR_TEST_OPTIONAL=true` when running
the command. Those tests are slower and depend on the external versions.

Help pages and `NAMESPACE` are maintained explicitly. Do not run an
unreviewed `roxygenise()` over them: several aliases and usage sections
are hand maintained. `Rscript tools/document.R` validates the help syntax.

Keep compiled objects, downloaded dependencies, datasets and experiment
outputs out of Git. Source code belongs in `R/` and `src/`; installed
examples in `inst/examples/`; regression tests in `tests/testthat/`.

When reporting a numerical issue, include a small reproducible input,
`packageVersion("egcar")`, `sessionInfo()`, the complete controls, and the
fit's status/convergence diagnostics. Open an issue at
<https://github.com/TannisthaM/EGCAR/issues>.

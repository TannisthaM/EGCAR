# Run explicitly in a fresh R session. Optional benchmarks are not installed here.
REPO <- "TannisthaM/EGCAR"
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")
# Installs the source package and compiles its C++ backend during installation.
remotes::install_github(REPO, dependencies = NA, upgrade = "never", build_vignettes = FALSE)
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.18")
print(benchmark_dependencies(), row.names = FALSE)

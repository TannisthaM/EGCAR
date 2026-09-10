# Run in the RStudio console, or edit REPO and source this file.
# Replace the placeholder with the repository you published through GitHub Desktop.
REPO <- "YOUR_GITHUB_USERNAME/egcar"
if (grepl("YOUR_GITHUB_USERNAME", REPO, fixed = TRUE)) {
  stop("Set REPO to your actual GitHub owner/repository before running this installer.")
}
packages <- c("remotes", "Rcpp", "RcppArmadillo", "future", "future.apply",
              "RGCCA", "PMA", "ggplot2", "RSpectra", "RhpcBLASctl")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# Installs the source package and compiles its C++ backend during installation.
remotes::install_github(REPO, dependencies = NA, upgrade = "never", build_vignettes = FALSE)
library(egcar)
stopifnot(packageVersion("egcar") >= "0.2.0")
print(benchmark_dependencies(), row.names = FALSE)

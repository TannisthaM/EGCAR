#!/usr/bin/env Rscript
# Help pages and NAMESPACE are shipped explicitly for this research revision.
# Validate them without overwriting hand-maintained aliases or usage sections.
for (f in list.files("man", pattern = "\\.Rd$", full.names = TRUE)) tools::parse_Rd(f)
cat("Help pages parsed. Use R CMD check for usage and namespace validation.\n")

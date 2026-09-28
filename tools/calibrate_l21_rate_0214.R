#!/usr/bin/env Rscript
# Run from the package root:
# Rscript tools/calibrate_l21_rate_0214.R pilot OUTPUT_DIRECTORY
# Rscript tools/calibrate_l21_rate_0214.R select OUTPUT_DIRECTORY
# Rscript tools/calibrate_l21_rate_0214.R confirmation OUTPUT_DIRECTORY
# Optional third argument: comma-separated scenario indices, for partitioned jobs.
# Set EGCAR_CALIBRATION_LIBRARY only to select an isolated installed package.
local_lib <- Sys.getenv("EGCAR_CALIBRATION_LIBRARY", "")
if (nzchar(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))
library(egcar)
source("inst/validation/l21_pilot_helpers.R")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L, args[[1L]] %in% c("pilot", "select", "confirmation"))
stage <- args[[1L]]
outdir <- args[[2L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
grid <- c(0.05, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75, 1.00)
if (stage == "select") {
  d <- l21_pilot_design("pilot")
  files <- file.path(outdir, sprintf("pilot_%03d.csv", d$scenario))
  if (!all(file.exists(files))) stop("Complete every prespecified pilot scenario before selecting.")
  x <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  stopifnot(nrow(x) == nrow(d) * length(grid), !anyDuplicated(x[c("scenario", "multiplier")]))
  summary <- do.call(rbind, lapply(grid, function(cg) {
    a <- x[x$multiplier == cg, ]
    good <- a$status == "ok" & a$converged & is.finite(a$normalized_score)
    eligible <- length(good) == nrow(d) && all(good)
    data.frame(multiplier = cg, scenarios = nrow(a), valid_converged = sum(good),
      eligible = eligible, mean_normalized_score = if (eligible) mean(a$normalized_score) else NA_real_,
      mean_selected_fraction = mean(a$selected_rows / (3 * a$p_per_block)))
  }))
  eligible <- which(summary$eligible)
  if (!length(eligible)) stop("No common multiplier is valid and converged in every pilot scenario.")
  # Prespecified rule: maximize mean held-out score / rank. Exact ties favor
  # the larger multiplier. No ground-truth errors enter this calculation.
  best <- eligible[order(-summary$mean_normalized_score[eligible],
                         -summary$multiplier[eligible])[[1L]]]
  summary$selected <- seq_len(nrow(summary)) == best
  selected <- summary$multiplier[[best]]
  frozen <- file.path(outdir, "selected_multiplier.txt")
  if (file.exists(frozen) && as.numeric(readLines(frozen)) != selected)
    stop("A different multiplier was already frozen; use a new study directory.")
  write.csv(x, file.path(outdir, "pilot_results.csv"), row.names = FALSE)
  write.csv(summary, file.path(outdir, "pilot_summary.csv"), row.names = FALSE)
  writeLines(format(selected, digits = 17), frozen)
  print(summary, row.names = FALSE)
  quit(save = "no")
}
d <- l21_pilot_design(stage)
write.csv(d, file.path(outdir, paste0(stage, "_design.csv")), row.names = FALSE)
if (stage == "confirmation") {
  frozen <- file.path(outdir, "selected_multiplier.txt")
  if (!file.exists(frozen)) stop("Freeze the pilot-selected multiplier before confirmation.")
  grid <- sort(unique(c(as.numeric(readLines(frozen)), 1)))
}
indices <- if (length(args) >= 3L) as.integer(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]]) else d$scenario
stopifnot(all(indices %in% d$scenario), !anyDuplicated(indices))
control <- egcar_control(max_iter = 2000L, abs_tol = 1e-5, rel_tol = 1e-4,
                         check_every = 5L, keep_full_C = FALSE, compact_state = TRUE)
capture.output(sessionInfo(), file = file.path(outdir, paste0(stage, "_session.txt")))
for (i in indices) {
  target <- file.path(outdir, sprintf("%s_%03d.csv", stage, i))
  if (file.exists(target)) next
  result <- l21_pilot_evaluate(d[i, ], grid, control)
  write.csv(result, paste0(target, ".tmp"), row.names = FALSE)
  if (!file.rename(paste0(target, ".tmp"), target)) stop("Could not commit scenario results.")
  invisible(gc(FALSE))
}

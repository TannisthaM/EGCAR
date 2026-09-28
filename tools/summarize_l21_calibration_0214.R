#!/usr/bin/env Rscript
# Run from the package root after all confirmation cases finish:
# Rscript tools/summarize_l21_calibration_0214.R OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
outdir <- args[[1L]]
source("inst/validation/l21_pilot_helpers.R")
design <- l21_pilot_design("confirmation")
selected <- as.numeric(readLines(file.path(outdir, "selected_multiplier.txt")))
files <- file.path(outdir, sprintf("confirmation_%03d.csv", design$scenario))
if (!all(file.exists(files))) stop("Confirmation is incomplete.")
x <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
stopifnot(nrow(x) == nrow(design) * 2L,
          !anyDuplicated(x[c("scenario", "multiplier")]),
          setequal(x$multiplier, c(selected, 1)))
summarize <- function(a) {
  good <- a$status == "ok" & a$converged & is.finite(a$normalized_score)
  distance <- a$subspace_sigma0[good] / sqrt(a$rank[good])
  data.frame(
    scenarios = nrow(a), valid_converged = sum(good),
    invalid_loading = sum(a$status == "invalid_loading"),
    not_converged = sum(a$status == "not_converged"),
    errors = sum(a$status == "error"),
    certified_zero = sum(a$zero_certified),
    admm_fits = sum(!a$zero_certified),
    mean_score_per_rank_valid = if (any(good)) mean(a$normalized_score[good]) else NA_real_,
    mean_selected_fraction_valid = if (any(good)) mean(a$selected_rows[good] / (3 * a$p_per_block[good])) else NA_real_,
    min_selected_fraction_valid = if (any(good)) min(a$selected_rows[good] / (3 * a$p_per_block[good])) else NA_real_,
    mean_normalized_subspace_error_valid = if (length(distance)) mean(distance) else NA_real_,
    max_normalized_subspace_error_valid = if (length(distance)) max(distance) else NA_real_)
}
by_multiplier <- do.call(rbind, lapply(sort(unique(x$multiplier)), function(cg)
  cbind(multiplier = cg, summarize(x[x$multiplier == cg, ]))))
by_dimension <- do.call(rbind, lapply(sort(unique(x$p_per_block)), function(p)
  do.call(rbind, lapply(sort(unique(x$multiplier)), function(cg)
    cbind(p_per_block = p, multiplier = cg,
          summarize(x[x$multiplier == cg & x$p_per_block == p, ]))))))
write.csv(x, file.path(outdir, "confirmation_results.csv"), row.names = FALSE)
write.csv(by_multiplier, file.path(outdir, "confirmation_summary.csv"), row.names = FALSE)
write.csv(by_dimension, file.path(outdir, "confirmation_by_dimension.csv"), row.names = FALSE)
print(by_multiplier, row.names = FALSE)
print(by_dimension, row.names = FALSE)

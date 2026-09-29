#!/usr/bin/env Rscript
# Compatibility launcher. From 0.2.16 the installed package is the single source
# of the current SGCA algorithm; the historical duplicate implementation is retired.
run_local_egcar_experiments <- function(
    output_dir="egcar_local_outputs", workers=1L, n_reps=1L,
    config=NULL, backend="cpp", smoke_test=identical(Sys.getenv("EGCAR_SMOKE_TEST","0"),"1")) {
  if (!requireNamespace("egcar",quietly=TRUE) || utils::packageVersion("egcar") < "0.2.16")
    stop("Install egcar >= 0.2.16 before using this launcher.")
  if (is.null(config)) config <- egcar::egcar_experiment_config()
  egcar::run_egcar_experiments(output_dir,workers,n_reps,config,backend,smoke_test)
}
if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_local_egcar_experiments(
    output_dir=if(length(args)>=1L) args[[1L]] else "egcar_local_outputs",
    workers=if(length(args)>=2L) as.integer(args[[2L]]) else 1L,
    n_reps=if(length(args)>=3L) as.integer(args[[3L]]) else 1L)
}

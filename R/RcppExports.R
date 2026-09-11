# Namespace-resolved registered native routine: safe in multisession workers.
egcar_native_solve <- function(context, state, controls, group, verbose = FALSE) {
  .Call(`_egcar_egcar_native_solve`, context, state, controls, group, verbose)
}

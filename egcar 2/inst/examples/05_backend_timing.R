library(egcar)
# Fixed iteration counts isolate arithmetic costs. These are NOT statistical
# performance results and do not show the speed of a full CV study.
sim <- egcar_simulate(n = 40, p_list = c(70, 60, 50), rank = 2, seed = 20)
prepared <- egcar_prepare(sim$views)
out <- list()
for (penalty in c("l11", "l21")) for (backend in c("cpp", "R", "reference")) {
  elapsed <- replicate(3, {
    gc(FALSE)
    unname(system.time(egcar_fit(prepared, 2, penalty, lambda = 0.02,
      control = egcar_control(backend = backend, max_iter = 100,
        abs_tol = 0, rel_tol = 0, adaptive_mu = FALSE, check_every = 5)))[["elapsed"]])
  })
  out[[length(out) + 1L]] <- data.frame(penalty = penalty, backend = backend,
                                       median_seconds = median(elapsed))
}
print(do.call(rbind, out), row.names = FALSE)

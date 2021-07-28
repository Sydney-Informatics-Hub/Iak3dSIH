library(here)

setUptests <- function(lmm.fit.selected) {
  #save main outputs to facilitate testing
  saveRDS(lmm.fit.selected, file = here("tests/run_results/lmm.fit.selected.rds"))
  
}
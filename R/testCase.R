library(here)

setUptests <- function() {
  #save main outputs to facilitate testing
  saveRDS(lmm.fit.selected, file = here("tests/run_results/lmm.fit.selected.rds"))
  saveRDS(vkVal, file = here("tests/run_results/vkVal.rds"))
  saveRDS(zkVal, file = here("tests/run_results/zkVal.rds"))
  
}
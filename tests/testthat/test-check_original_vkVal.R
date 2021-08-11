library(testthat)

# Model run intentionally doesnt save large files to memory as transitioning to
# R package dev with limited memory usage. For these tests to run, you must perform 
# the saving manually.

# If Splinedata <- SplineIAK() run, then:
#saveRDS(Splinedata$xkvkVal$vkVal,file = here::here("tests/run_results/vkVal.rds"))
#saveRDS(Splinedata$xkvkVal$zkVal,file = here::here("tests/run_results/zkVal.rds"))
#saveRDS(Splinedata$lmm.fit.selected,file = here::here("tests/run_results/lmm.fit.selected.rds"))

# If Cubistdata <- CubistIAK() run then:
#saveRDS(Cubistdata$xkvkVal$vkVal,file = here::here("tests/run_results/vkVal.rds"))
#saveRDS(Cubistdata$xkvkVal$zkVal,file = here::here("tests/run_results/zkVal.rds"))
#saveRDS(Cubistdata$lmm.fit.selected,file = here::here("tests/run_results/lmm.fit.selected.rds"))


test_that("testing Spline model when run when fitCubistModelNow flag is FALSE", {
  
  expect_equal(readRDS(here::here("tests/original_results/vkVal.rds")),
               readRDS(here::here("tests/run_results/vkVal.rds")))
  
  expect_equal(readRDS(here::here("tests/original_results/zkVal.rds")),
               readRDS(here::here("tests/run_results/zkVal.rds")))
  
  
  expect_equal(readRDS(here::here("tests/original_results/lmm.fit.selected.rds")),
               readRDS(here::here("tests/run_results/lmm.fit.selected.rds")))
  
})

#only run this is fitCubistModelNow is TRUE in run
test_that("testing Cubist model when run when fitCubistModelNow flag is TRUE", {
  
  expect_equal(readRDS(here::here("tests/original_results_cubist/vkVal.rds")),
               readRDS(here::here("tests/run_results/vkVal.rds")))
  
  expect_equal(readRDS(here::here("tests/original_results_cubist/zkVal.rds")),
               readRDS(here::here("tests/run_results/zkVal.rds")))
  
  expect_equal(readRDS(here::here("tests/original_results_cubist/lmm.fit.selected.rds")),
               readRDS(here::here("tests/run_results/lmm.fit.selected.rds")))
})



             
library(testthat)

#only run this is fitCubistModelNow is FALSE in run
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



             
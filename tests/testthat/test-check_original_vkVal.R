library(testthat)

testData <- function(strobjrun,strobktest) {
  objThatRan <- readRDS(strobjrun)
  objTest <- readRDS(strobktest)
  identical(objThatRan,objTest)
}

test_that("Main RData files generated as output unchanged - Spline Model", {
  
  expect_equal(2 * 2, 4)
  expect_equal(testData(here::here("tests/original_results/vkVal.rds"),here::here("tests/run_results/vkVal.rds")),TRUE)
  expect_equal(testData(here::here("tests/original_results/zkVal.rds"),here::here("tests/run_results/zkVal.rds")),TRUE)
  expect_equal(testData(here::here("tests/original_results/lmm.fit.selected.rds"),here::here("tests/run_results/lmm.fit.selected.rds")),TRUE)
})

test_that("Main RData files generated as output unchanged Cubist Model", {
  expect_equal(testData(here::here("tests/original_results_cubist/vkVal.rds"),here::here("tests/run_results/vkVal.rds")),TRUE)
  expect_equal(testData(here::here("tests/original_results_cubist/zkVal.rds"),here::here("tests/run_results/zkVal.rds")),TRUE)
  expect_equal(testData(here::here("tests/original_results_cubist/lmm.fit.selected.rds"),here::here("tests/run_results/lmm.fit.selected.rds")),TRUE)
})

test_that("zkVal Cubist Model", {
  expect_equal(testData(here::here("tests/original_results_cubist/zkVal.rds"),here::here("tests/run_results/zkVal.rds")),TRUE)
})

test_that("vkVal Cubist Model", {
  expect_equal(testData(here::here("tests/original_results_cubist/vkVal.rds"),here::here("tests/run_results/vkVal.rds")),TRUE)
})

test_that("lmm.fit.selected.rdsCubist Model", {
  expect_equal(testData(here::here("tests/original_results_cubist/lmm.fit.selected.rds"),here::here("tests/run_results/lmm.fit.selected.rds")),TRUE)
})



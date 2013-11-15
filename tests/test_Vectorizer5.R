context("Vector checks")

test_that("SufficientStats works", {
  expect_that(SufficientStats(kmer.list = list(c(0,1)), j = 3),
    equals(0.5))
})

test_that("SufficientStats works", {
  expect_that(SufficientStats(kmer.list = list(c(0,1)), j = 3),
              equals(0.5))
})

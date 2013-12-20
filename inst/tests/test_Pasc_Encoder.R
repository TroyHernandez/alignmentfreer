# test_Pasc_Encoder.R
context("Pasc_Encoder checks")

test_that("Pasc_Encoder works", {
  
  target.mat <- read.csv("Alloherpesviridae.csv", header = TRUE, row.names = 1)
  colnames(target.mat) <- rownames(target.mat)
  target.mat <- as.matrix(target.mat)
  test.mat <- "Alloherpesviridae.txt"
  expect_that(print(pasc(test.mat)),
              equals(target.mat))
})


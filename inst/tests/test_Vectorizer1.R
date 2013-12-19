# test_Vectorizer1.R
context("Vectorizer1 checks")

test_that("Vectorizer works", {
  expect_that(as.numeric(Vectorizer(dna.seq = "GATTACA", kmer = 1)),
              equals(c(7, 3 / 7, 1 / 7, 1 / 7, 2 / 7,
                       2 / 3, 6 / 7, 1 / 7, 1 / 2,
                       var(c(2, 5, 7) / 7), 0, 0, var(c(3, 4) / 7) )))
})

test_that("Vectorizer works", {
  expect_that(as.numeric(Vectorizer(dna.seq = "NA", kmer = 1)),
              equals(c(2, 1.25 / 2, .25 / 2, .25 / 2, .25 / 2,
                       1.8, 1, 1, 1, 0.5, 0, 0, 0)))
})


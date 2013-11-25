# test_Vectorizer3.R

test_that("KmerGenerator works", {
  expect_that(KmerGenerator(kmer = 1),
              equals(c("A", "C", "G", "T")))
})

test_that("CalcRegLetters works", {
  expect_that(CalcRegLetters(kmer.seq = c("G", "A", "T", "T", "A", "C", "A"),
                             tbl = c(A = 3, C = 1, G = 1, T = 2),
                             statistic = 3, kmers = c("A", "C", "G", "T"),
                             method = "Sufficient"),
              equals(c(3 / 7, 1 / 7, 1 / 7, 2 / 7,
                       2 / 3, 6 / 7, 1 / 7, 1 / 2,
                       var(c(2, 5, 7) / 7), 0, 0, var(c(3, 4) / 7) )))
})

test_that("CalcRegLetters works", {
  expect_that(CalcRegLetters(kmer.seq = c("G", "A", "T", "T", "A", "C", "A"),
                             tbl = c(A = 3, C = 1, G = 1, T = 2),
                             statistic = 3, kmers = c("A", "C", "G", "T"),
                             method = "Sufficient"),
              equals(c(3 / 7, 1 / 7, 1 / 7, 2 / 7,
                       2 / 3, 6 / 7, 1 / 7, 1 / 2,
                       var(c(2, 5, 7) / 7), 0, 0, var(c(3, 4) / 7) )))
})

test_that("CalcAmbigLetters works", {
  expect_that(CalcAmbigLetters(kmer.seq = c("N", "A"),
                             statistic = 3, kmers = c("A", "C", "G", "T"),
                             method = "Sufficient"),
              equals(c(5 / 8, 1 / 8, 1 / 8, 1 / 8,
                       9 / 5, 1, 1, 1,
                       .5, 0, 0, 0)))
})

test_that("CalcEmptyLetters works", {
  expect_that(CalcEmptyLetters(tbl, statistic, kmers, method = "Sufficient"),
              equals(c(5 / 8, 1 / 8, 1 / 8, 1 / 8,
                       9 / 5, 1, 1, 1,
                       .5, 0, 0, 0)))
})

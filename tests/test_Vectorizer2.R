# test_Vectorizer2.R

test_that("SanityCheck works", {
  expect_that(SanityCheck(dna.seq = "gattaca"),
              equals(NULL))
})

test_that("UpperCaser works", {
  expect_that(UpperCaser(dna.seq = "gattaca"),
              equals("GATTACA"))
})

test_that("KmerColNames works", {
  expect_that(KmerColNames(kmer = 2, statistic = 1, concatenate = F),
              equals(paste("n", c(t(outer(c("A", "C", "G", "T"),
                                          c("A", "C", "G", "T"),
                                          FUN = "paste", sep = ""))),
                           sep = "")))
})

test_that("CalculateVec works", {
  expect_that(CalculateVec(c("G", "A", "T", "T", "A", "C", "A"),
                           statistic = 3, kmer = 1, method = "Sufficient"),
              equals(list(vector = c(3 / 7, 1 / 7, 1 / 7, 2 / 7,
                       2 / 3, 6 / 7, 1 / 7, 1 / 2,
                       var(c(2, 5, 7) / 7), 0, 0, var(c(3, 4) / 7) ),
                          length = 7)))
})

ColumnFinder(k = 2, statistic = 2, length.vec = 1, concatenate = F)
test_that("ColumnFinder works", {
  expect_that(ColumnFinder(k = 2, statistic = 2,
                           length.vec = 1, concatenate = T),
              equals(10:41))
})

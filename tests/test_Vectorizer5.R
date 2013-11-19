context("Vector checks")

test_that("SufficientStats works", {
  expect_that(SufficientStats(kmer.list = list(c(0,1)), j = 3),
    equals(0.5))
})

test_that("AllocateIrregKmers works", {
  expect_that(AllocateIrregKmers(irreg.name = "N",
                                 kmers = c("A", "C", "G", "T"),
                                 reg.tbl = rep(0,4), tbl = table("N"),
                                 kmer.seq = "N",
                                 kmer.list = list(A = 0, C = 0, G = 0, T = 0),
                                 kmer.wt.list = list(A = 0, C = 0,
                                                     G = 0, T = 0)),
              equals(list(kmer.list = list(A = c(0, 1), C = c(0, 1),
                                           G = c(0, 1), T = c(0, 1)),
                          kmer.wt.list = list(A = c(0, .25), C = c(0, .25),
                                           G = c(0, .25), T = c(0, .25)),
                          reg.tbl = rep(.25, 4))))
})

test_that("IrregSufficientStats works", {
  expect_that(IrregSufficientStats(kmer.list = list(c(0, 1)),
                                   kmer.wt.list = list(c(.25, 1)),
                                   statistic = 3),
              equals(0.5))
})

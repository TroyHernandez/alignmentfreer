context("Vector checks")

test_that("CalcNumLocIrreg works", {
  expect_that(CalcNumLocIrreg("D"),
    equals(list(num.perms = 3, extra.letters.vec = 2)))
})

test_that("CalcNumLocIrreg works", {
  expect_that(CalcNumLocIrreg("N"),
              equals(list(num.perms = 4, extra.letters.vec = 1)))
})

test_that("CalcPermutationMat works", {
  expect_that(CalcPermutationMat("RY", c(2, 2)),
              equals(matrix(c("G", "G", "A", "A", "T", "C", "T", "C"),
                            nrow = 4)))
})

test_that("AddToTblKmer works", {
  expect_that(AddToTblKmer(perm.mat = CalcPermutationMat("N",4),
                           kmers = c("A", "C", "G", "T"), reg.tbl = rep(0,4),
                           tbl = table("N"), kmer.seq = "N",
                           kmer.list = list(A = 0, C = 0, G = 0, T = 0),
                           kmer.wt.list = list(A = 0, C = 0, G = 0, T = 0),
                           irreg.name = "N"),
              equals(list(kmer.list = list(A = c(0, 1), C = c(0, 1),
                                           G = c(0, 1), T = c(0, 1)),
                          kmer.wt.list = list(A = c(0, .25), C = c(0, .25),
                                           G = c(0, .25), T = c(0, .25)),
                          reg.tbl = rep(0.25,4))))
})

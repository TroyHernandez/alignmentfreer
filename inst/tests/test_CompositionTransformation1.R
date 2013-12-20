# test_CompositionTransformation1.R
context("CompTrans checks")

test_that("CompositionTransformations works", {
  vec <- Vectorizer(dna.seq = "GATTACA", kmer = 2, statistic = 1)
  
  kGattaca1Mer <- c(3 / 7, 1 / 7, 1 / 7, 2 / 7)
  kExp2Mer <- outer(kGattaca1Mer, kGattaca1Mer)
  
  test.value <- CompTransform(vec, kmer = 2, statistic = 1)
  target.value <- c(7, 3 / 7, 1 / 7, 1 / 7, 2 / 7,
                    (0 / 6 - kExp2Mer[1, 1]) / sqrt(kExp2Mer[1, 1]),
                    (1 / 6 - kExp2Mer[1, 2]) / sqrt(kExp2Mer[1, 2]),
                    (0 / 6 - kExp2Mer[1, 3]) / sqrt(kExp2Mer[1, 3]),
                    (1 / 6 - kExp2Mer[1, 4]) / sqrt(kExp2Mer[1, 4]),
                    (1 / 6 - kExp2Mer[2, 1]) / sqrt(kExp2Mer[2, 1]),
                    (0 / 6 - kExp2Mer[2, 2]) / sqrt(kExp2Mer[2, 2]),
                    (0 / 6 - kExp2Mer[2, 3]) / sqrt(kExp2Mer[2, 3]),
                    (0 / 6 - kExp2Mer[2, 4]) / sqrt(kExp2Mer[2, 4]),
                    (1 / 6 - kExp2Mer[3, 1]) / sqrt(kExp2Mer[3, 1]),
                    (0 / 6 - kExp2Mer[3, 2]) / sqrt(kExp2Mer[3, 2]),
                    (0 / 6 - kExp2Mer[3, 3]) / sqrt(kExp2Mer[3, 3]),
                    (0 / 6 - kExp2Mer[3, 4]) / sqrt(kExp2Mer[3, 4]),
                    (1 / 6 - kExp2Mer[4, 1]) / sqrt(kExp2Mer[4, 1]),
                    (0 / 6 - kExp2Mer[4, 2]) / sqrt(kExp2Mer[4, 2]),
                    (0 / 6 - kExp2Mer[4, 3]) / sqrt(kExp2Mer[4, 3]),
                    (1 / 6 - kExp2Mer[4, 4]) / sqrt(kExp2Mer[4, 4]))
  
  names(target.value) <- names(test.value)
  
  expect_that(as.numeric(test.value), equals(target.value))
})

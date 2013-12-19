#test_Vectorizer4.R
context("Vectorizer4 checks")

test_that(".CalcMeanKmerList works", {
  expect_that(.CalcMeanKmerList(kmer.seq = c("G", "A", "T", "T", "A", "C", "A"),
                               tbl = table(c("G", rep("A", 3),
                                             rep("T", 2), "C")),
                               statistic = 3,
                               kmers = c("A", "C", "G", "T")),
              equals(list(mean = c(14 / 21, 6 / 7, 1 / 7, 7 / 14),
                          list = list(c(2 / 7, 5 / 7, 7 / 7), c(6 / 7),
                                      c(1 / 7), c(3 / 7, 4 / 7)))))
})

test_that(".CalcDescriptiveStats works", {
  expect_that(.CalcDescriptiveStats(ans = matrix(rep(0, 12), nrow = 3),
                                   kmer.list = list(c(2 / 7, 5 / 7, 7 / 7),
                                                    c(6 / 7), c(1 / 7),
                                                    c(3 / 7, 4 / 7)),
                               statistic = 3,
                               method = "Sufficient")[3, ],
              equals(c(var(c(2 / 7, 5 / 7, 7 / 7)), NA,
                       NA, var( c(3 / 7, 4 / 7)))))
})

test_that(".CalcAmbigKmerList works", {
  expect_that(.CalcAmbigKmerList(kmer.seq = c("N", "A"),
                                tbl = c(A = 1, N = 1),
                                ambig.names = c("A", "N"),
                                kmers = c("A", "C", "G", "T"),
                                ans = matrix(0, nrow = 3, ncol = 4)),
              equals(list(kmer.list = list(A = c(2, 1), C = 1, G = 1, T = 1),
                          kmer.wt.list = list(A = c(1, .25), C = .25,
                                              G = .25, T = .25),
                          tbl = c(A = 1.25, C = .25, G = .25, T = .25),
                          ans = matrix(c(1.25, 0, 0, .25, 0, 0,
                                         .25, 0, 0, .25, 0, 0), nrow = 3))))
})

test_that(".CalcAmbigDescriptiveStats works", {
  expect_that(.CalcAmbigDescriptiveStats(ans = matrix(c(1.25, 0, 0, .25, 0, 0,
                                                       .25, 0, 0, .25, 0, 0),
                                                     nrow = 3),
                                        kmer.list = list(A = c(2, 1), C = 1,
                                                         G = 1, T = 1),
                                        kmer.wt.list = list(A = c(1, .25),
                                                            C = .25, G = .25,
                                                            T = .25),
                                        statistic = 3, method = "Sufficient"),
              equals(matrix(c(1.25, 1.8, 0.5, .25, 1, NaN,
                              .25, 1, NaN, .25, 1, NaN), nrow = 3)))
})


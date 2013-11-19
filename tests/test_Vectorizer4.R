#test_Vectorizer4.R

test_that("CalcMeanKmerList works", {
  expect_that(CalcMeanKmerList(kmer.seq = c("G", "A", "T", "T", "A", "C", "A"),
                               tbl = table(c("G", rep("A", 3),
                                             rep("T", 2), "C")),
                               statistic = 3,
                               kmers = c("A", "C", "G", "T")),
              equals(list(mean = c(14 / 21, 6 / 7, 1 / 7, 7 / 14),
                          list = list(c(2 / 7, 5 / 7, 7 / 7), c(6 / 7),
                                      c(1 / 7), c(3 / 7, 4 / 7)))))
})

test_that("CalcDescriptiveStats works", {
  expect_that(CalcDescriptiveStats(ans = matrix(rep(0, 12), nrow = 3),
                                   kmer.list = list(c(2 / 7, 5 / 7, 7 / 7),
                                                    c(6 / 7), c(1 / 7),
                                                    c(3 / 7, 4 / 7)),
                               statistic = 3,
                               method = "Sufficient")[3, ],
              equals(c(var(c(2 / 7, 5 / 7, 7 / 7)), NA,
                       NA, var( c(3 / 7, 4 / 7)))))
})


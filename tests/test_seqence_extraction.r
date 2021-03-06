context("Sequence extraction checks")

test_that("gbk extraction works", {
  expect_that(as.character(GetSeqGbk("data/abalone.gbk")),
              equals(as.character(read.table("data/abalone_string.txt",
                                             stringsAsFactors = FALSE,
                                             header = TRUE))))
})

test_that("fasta extraction works", {
  expect_that(as.character(GetSeqFasta("data/abalone.fasta", lower = TRUE)),
              equals(as.character(read.table("data/abalone_string.txt",
                                             stringsAsFactors = FALSE,
                                             header = TRUE))))
})

test_that("gbk extraction equivalent fasta", {
  fasta.test <- GetSeqFasta("data/abalone.fasta")
  gbk.test <- GetSeqGbk("data/abalone.gbk", upper = TRUE)
  expect_that(as.character(fasta.test),
              equals(as.character(gbk.test)))
  expect_that(attributes(fasta.test)$accession,
              equals(attributes(gbk.test)$accession))
  expect_that(attributes(fasta.test)$gi,
              equals(attributes(gbk.test)$gi))
  expect_that(attributes(fasta.test)$bp,
              equals(attributes(gbk.test)$bp))
})

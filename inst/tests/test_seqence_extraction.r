context("Sequence extraction checks")

test_that("gbk extraction works", {
  expect_that(as.character(gbk("abalone.gbk")),
              equals(read.table("abalone_string.txt", stringsAsFactors = FALSE,
                                header = TRUE)[[1]]))
})

test_that("fasta extraction works", {
  expect_that(as.character(fasta("abalone.fasta", lower = T)),
              equals(read.table("abalone_string.txt",
                                stringsAsFactors = FALSE, header = TRUE)[[1]]))
})

test_that("gbk extraction equivalent fasta", {
  fasta.test <- fasta("abalone.fasta")
  gbk.test <- gbk("abalone.gbk", upper = T)
  expect_that(as.character(fasta.test),
              equals(as.character(gbk.test)))
  expect_that(attributes(fasta.test)$accession,
              equals(attributes(gbk.test)$accession))
  expect_that(attributes(fasta.test)$gi,
              equals(attributes(gbk.test)$gi))
  expect_that(attributes(fasta.test)$bp,
              equals(attributes(gbk.test)$bp))
})

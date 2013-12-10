# Vectorizer3.R

#' Creates a vector of kmers for a given k value.
#'
#' This function creates a vector of length 4^k consisting of the kmers for a
#' given value of k. 
#'
#' @param kmer an integer that defines the kmer length to be calculated 
#' @keywords vectorizer
#' @export
#' 
KmerGenerator <- function(kmer) {
  
  dna <- c("A","C","G","T")
  
  #kmercolnames is final output
  kmer.list <- dna
  
  if (kmer > 1) {
    for (j in 2:kmer) {
      kmer.list <- c(t(outer(kmer.list, dna, FUN = "paste", sep = "")))
    }
  }
  
  kmer.list
}

#######################################################
# Appends statistic name to kmer 
.StatAppender <- function(kmers, statistic) {
  for (i in 1:statistic) {
    if (i == 1) {
      col.names <- paste("n", kmers, sep = "")
    }
    if (i == 2) {
      col.names <- c(col.names,paste("mu", kmers, sep = ""))
    }
    if (i == 3) {
      col.names <- c(col.names,paste("var", kmers, sep = ""))
    }
    if (i == 4) {
      col.names <- c(col.names,paste("skew", kmers, sep = ""))
    }
    if (i == 5) {
      col.names <- c(col.names,paste("kur", kmers, sep = ""))
    }
  }
  col.names
}

##############################################################
.CalcRegLetters <- function(kmer.seq, tbl, statistic, kmers,
                           method = "Sufficient") {
  ans <- matrix(0, nrow = statistic, ncol = length(kmers))
  ans[1, ] <- tbl
  
  if (statistic > 1) {
    mean.kmer.list <- .CalcMeanKmerList(kmer.seq, tbl, statistic, kmers)
    ans[2, ] <- mean.kmer.list$mean
    kmer.list <- mean.kmer.list$list
  }
  ans <- .CalcDescriptiveStats(ans, kmer.list, statistic, method)
  
  if (statistic > 1) {
    ans <- .CorrectZeroCases(ans, statistic)
    ans <- .CorrectSingletonCases(ans, statistic)
  }
  
  # Changes counts to frequencies
  ans[1, ] <- ans[1, ] / sum(ans[1, ])
  
  ans <- c(t(ans))
}



################################################################

.CalcAmbigLetters <- function (kmer.seq, statistic, kmers,
                              method = "Sufficient") {
  ans <- matrix(0, nrow = statistic, ncol = length(kmers))
  tbl <- table(kmer.seq)
  ambig.names <- names(tbl)
  
  ambig.kmer.list <- .CalcAmbigKmerList(kmer.seq, tbl,
                                           ambig.names, kmers, ans)
  tbl <- ambig.kmer.list$tbl
  kmer.list <- ambig.kmer.list$kmer.list
  kmer.wt.list <- ambig.kmer.list$kmer.wt.list
  ans <- ambig.kmer.list$ans

  ans <- .CalcAmbigDescriptiveStats(ans, kmer.list, kmer.wt.list,
                                   statistic, method)
  
  if (statistic > 1) {
    ans <- .CorrectZeroCases(ans, statistic)
    ans <- .CorrectSingletonCases(ans, statistic)
  }
  
  ans[1, ] <- ans[1, ] / sum(ans[1, ])
  
  ans <- c(t(ans))
  ans
}

##############################################################
.CalcEmptyLetters=function(tbl, statistic, kmers, method = "Sufficient"){
  ans <- matrix(0, nrow = statistic, ncol = length(kmers))
  ans <- .CorrectZeroCases(ans, statistic, method)
  ans <- c(t(ans))
  ans
}

# Vectorizer.R


library(e1071)
source("Vectorizer2.R")
source("Vectorizer3.R")
source("Vectorizer4.R")
source("Vectorizer5.R")
source("Vectorizer6.R")
#######################################################
#Main function for calculating natural vectors
Vectorizer <- function(dna.seq, kmer = 1, statistic = 3,
                       composition = TRUE, concatenate = TRUE,
                       confirm.dna.seq = FALSE, verbose = 0) {
  
  if (confirm.dna.seq == TRUE) {
    ConfirmDnaSeq(dna.seq)
  }
  dna.seq <- UpperCaser(dna.seq)
  range <- .Ranger(kmer, concatenate)
  col.names <- .KmerColNames(kmer, statistic, concatenate)
  vec <- matrix(0, nrow = 1, ncol = length(col.names) + 1)
  colnames(vec) <- c("length", col.names)
  vec[1] <- nchar(dna.seq)
  
  for (k in range) {
    if(verbose == 1) {
      cat(".")      
    }
    
    if(verbose == 2) {
      cat("Calculating ", k, "-mer\n", sep = "")      
    }
#     Test string
#     dna.seq <- "NA"
#     dna.seq <- ""
#     dna.seq <- "GATTACA"
#     dna.seq <- paste(sample(c("A","C","G","T","B","N"),12,replace=T),
#                      collapse="")
    
    num.kmer <- nchar(dna.seq) - k + 1
    if (num.kmer > 0) {
      kmer.seq <- substring(dna.seq, 1:num.kmer, k:nchar(dna.seq))
    } else {
      kmer.seq <- ""
    }
    temp.vec <- CalculateVec(kmer.seq, statistic,
                             kmer = k, method = "Sufficient")
    
    temp.cols <- .ColumnFinder(k, statistic, length(vec), concatenate)
    vec[1, temp.cols] <- temp.vec$vector
  }
  vec
}

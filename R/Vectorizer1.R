# Vectorizer.R
#=========================================================

library(e1071)
source("Vectorizer2.R")
source("Vectorizer3.R")
source("Vectorizer4.R")
source("Vectorizer5.R")
source("Vectorizer6.R")
#=========================================================

#' Calculate a generalized vector from a DNA sequence.
#'
#' This function takes a DNA sequence in character string format; e.g. "GATTACA" 
#' and calculates alignment free statistics.
#'
#' @param dna.seq an upper or lower-case character vector consisting of A, C, T,
#'   G, or any of the IUPAC ambiguous nucleotide letters
#' @param kmer an integer that defines the kmer length to be calculated 
# ' @param composition logical indicating whether kmer frequencies should be
# '   returned if the composition vector should be returned
#' @param concatenate logical indicating whether the statistics for the kmer
#'   should be returned or if the 1st through kth kmer statistics should be 
#'   returned
#' @param confirm.dna.seq logical indicating that the \code{dna.seq} should be
#'   confirmed as a proper dna sequence; this will add to the run-time
#' @param verbose integer taking the values 0, 1, or 2 for varying levels of
#'   verbosity
#' @keywords vectorizer
#' @export
#' @examples
#' 
#' 
#' # K-mer with k = 2
#' Vectorizer(GetSeqGbk("data/abalone.gbk", upper = T),
#'                       kmer = 2, statistic = 1, concatenate = F)
#' #      length       nAA        nAC        nAG        nAT       nCA        nCC
#' # [1,]  34952 0.1059197 0.05338903 0.06194386 0.08572001 0.0779663 0.04094303
#' #            nCG        nCT        nGA        nGC        nGG        nGT        nTA
#' # [1,] 0.0219164 0.06283082 0.05267374 0.05130039 0.03705187 0.05052788 0.07044148
#' #           nTC        nTG        nTT
#' [1,] 0.05799548 0.07064176 0.09873823
#' 
#' # The Natural Vector
#' abalone <- Vectorizer(GetSeqGbk("data/abalone.gbk", upper = T),
#'                       kmer = 1, statistic = 3)
#' abalone.natural.vector <- abalone[1]*abalone[2:13]
#' abalone.natural.vector
#' # [1] 10730.000  7118.000  6695.000 10409.000 17319.117 17868.755 17451.659
#' # [8] 17386.478  2946.577  2794.014  2914.679  2954.686#' 
#' 
#' 
#' 
Vectorizer <- function(dna.seq, kmer = 3, statistic = 3,
#                        composition = TRUE,
                       concatenate = TRUE,
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
    
    if(verbose >= 2) {
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
    temp.vec <- .CalculateVec(kmer.seq, statistic,
                             kmer = k, method = "Sufficient")
    
    temp.cols <- .ColumnFinder(k, statistic, length(vec), concatenate)
    vec[1, temp.cols] <- temp.vec$vector
  }
  vec
}

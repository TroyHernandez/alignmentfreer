# CompositionTransformation.R
#=========================================================

#' Calculate the composition vector from a generalized vector sequence.
#'
#' This function takes a generalized vector sequence, as obtained from
#' \code{Vectorizer}, and calculates the composition vector.
#'
#' @param vec a generalized vector sequence as derived from \code{Vectorizer}
#' with kmer greater than 1.
#' @param kmer an integer that defines the kmer length of vec.
#' @param statistic an integer that defines the statistics present in vec.
#' @keywords vectorizer
#' @export
#' 
#' 
CompTransform <- function(vec, kmer = 3, statistic = 3) {
    
  if (kmer < 2) {
    stop("This function only works for kmers with k > 1!")
  }
  
  #Work backwards through kmers
  for (i in kmer:2) {
    for (j in 1:(4 ^ i)) {
      #target kmer with "n" in the front
      nkmer <- colnames(vec)[.ColumnFinder(k = kmer, statistic = statistic,
                                       length.vec = length(vec),
                                        concatenate = T)][j]
      #target kmer without "n" in the front
      temp.kmer <- substring(nkmer, 2, i + 1)
      #target left sub-kmer with "n" in the front
      n.l.kmer <- paste("n", substring(temp.kmer, 1, i - 1), sep = "")
      #target right sub-kmer with "n" in the front
      n.r.kmer <- paste("n", substring(temp.kmer, 2, i), sep = "")
      #target center sub-kmer with "n" in the front
      if (i > 2) {
        n.w.kmer <- paste("n", substring(temp.kmer, 2, i - 1), sep = "")
        n.w.kmer.vec <- vec[, n.w.kmer]
      }else{
        n.w.kmer.vec <- 1
      }
      #Create q vector
      qvec <- vec[, n.l.kmer] * vec[, n.r.kmer] / n.w.kmer.vec
      #Replace with composition vector
      vec[, nkmer] <- (vec[, nkmer] - qvec) / sqrt(qvec)
      
      #Check for qvec==0
      badq <- which(is.na(vec[, nkmer]) | abs(vec[, nkmer]) == Inf)
      if (length(badq) > 0) {
        vec[badq, nkmer] <- 0
      }
    }
  }
  vec
}

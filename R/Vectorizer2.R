# Vectorizer2.R
#=========================================================

#' Confirms character is a valid DNA sequence.
#'
#' This function checks that a DNA sequence in character string format,
#' e.g. "GATTACA", and confirms it consists only of upper or lower-case
#' character vector consisting of A, C, T, G, or any of the IUPAC ambiguous
#' nucleotide letters.
#'
#' @param dna.seq an upper or lower-case character vector consisting of A, C, T,
#'   G, or any of the IUPAC ambiguous nucleotide letters
#' @keywords vectorizer
#' @export
#' 
ConfirmDnaSeq <- function(dna.seq){
  dna.seq.letters <- strsplit(dna.seq, "")[[1]]
  k.AmbigLettersBig <- c("R", "Y", "M", "K", "S", "W", "B", "D", "H", "V", "N")
  k.AmbigLettersSml <- c("r", "y", "m", "k", "s", "w", "b", "d", "h", "v", "n")
  k.NucleotidesBig <- c("A", "C", "G", "T")
  k.NucleotidesSml <- c("a", "c", "g", "t")
  k.All.Letters <- c(k.AmbigLettersBig, k.AmbigLettersSml,
                     k.NucleotidesBig, k.NucleotidesSml)
  alphabet.len <- length(union(dna.seq.letters, k.All.Letters))
  if (alphabet.len > 30) {
    stop("Error: Nucleotide not recognized.\n
        Only A, C, G, T, or ambiguous letters allowed.")
  }
  alphabet.len.big <- length(union(dna.seq.letters,
                             c(k.AmbigLettersBig, k.NucleotidesBig)))
  alphabet.len.sml <- length(union(dna.seq.letters,
                             c(k.AmbigLettersSml, k.NucleotidesSml)))
  
  if (alphabet.len.big != 15 & alphabet.len.sml != 15) {
    stop("Error: Mixture of upper and lower case letters used.
        This is very atypical.  Check your data.")
  }
}

#=========================================================

#' Converts a character string to upper-case letters.
#'
#' This function converts a character string of letters to an upper-case
#' character string of letters.
#'
#' @param dna.seq an upper or lower-case character vector consisting of A, C, T,
#'   G, or any of the IUPAC ambiguous nucleotide letters
#' @keywords vectorizer
#' @export
#' 
UpperCaser <- function(dna.seq){
  if (sum(substring(dna.seq, 1, 1) == letters) > 0){
    ltrs=list(LETTERS)[[1]]
    names(ltrs)=letters
    dna.seq <- paste(ltrs[strsplit(dna.seq, "")[[1]]], collapse = "")
  }
  dna.seq
}

#=========================================================

# Determines range for (un)concatenated vectors.
.Ranger <- function(kmer, concatenate)  {
  if (concatenate) {
    range <- 1:kmer
  } else {
    range <- kmer
  }
}

#=========================================================

.KmerColNames <- function(kmer = 3, statistic = 3, concatenate = TRUE) {
  
  #kmercolnames is final output
  kmer.colnames <- c()
  
  range <- .Ranger(kmer, concatenate)
  
  for (j in range) {
    #append new round to total
    kmer.colnames <- c(kmer.colnames, .StatAppender(KmerGenerator(j), statistic))
  }
  
  kmer.colnames
}

.CalculateVec <- function(kmer.seq, statistic = 3,
                         kmer = 3, method = "Sufficient") {
  
  seq.length <- length(kmer.seq)
  ans <- rep(0, 4 ^ kmer * statistic)  
  kmers <- KmerGenerator(kmer)
  
  #This is used to ensure there are only A, C, G, T.
  tbl <- table(c(kmer.seq, kmers)) - 1
  
  #There aren't/are other letters
  if (length(tbl) == length(kmers)) {
    if (sum(names(tbl) != kmers) == 0) {
      ans <- .CalcRegLetters(kmer.seq, tbl, statistic, kmers, method)
    } else {
      ans <- .CalcAmbigLetters(kmer.seq, statistic, kmers, method)
    }
  } else {
    if (sum(kmer.seq == "") > 0) {
      ans <- .CalcEmptyLetters(tbl, statistic, kmers)
    } else {
      ans <- .CalcAmbigLetters(kmer.seq, statistic, kmers, method)
    }
  }
  
  list(vector = ans, length = seq.length);
}

#=========================================================

.ColumnFinder <- function(k, statistic, length.vec, concatenate) {
  # concatenate == F just removes the length column.
  if (!concatenate) {
    temp.cols <- 2:length.vec
  } else {
    if (k == 1) {
      temp.cols <- 1:(statistic * 4 ^ k)
      temp.cols <- temp.cols + 1
    } else {
      temp <- statistic * sum(4 ^ (1:(k - 1)))
      temp.cols <- (temp + 1):(temp + statistic * 4 ^ k)
      temp.cols <- temp.cols + 1
    }
  }
  temp.cols
}

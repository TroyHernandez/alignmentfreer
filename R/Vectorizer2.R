# Vectorizer2.R

UpperCaser <- function(dna.seq){
  if (sum(substring(dna.seq, 1, 1) == letters) > 0){
    ltrs=list(LETTERS)[[1]]
    names(ltrs)=letters
    dna.seq <- paste(ltrs[strsplit(dna.seq, "")[[1]]], collapse = "")
  }
  dna.seq
}
#######################################################
# Determines range for (un)concatenated vectors.
Ranger <- function(kmer, concatenate)  {
  if (concatenate) {
    range <- 1:kmer
  } else {
    range <- kmer
  }
}

#######################################################
KmerColNames <- function(kmer = 3, statistic = 3, concatenate = TRUE) {
  
  #kmercolnames is final output
  kmer.colnames <- c()
  
  range <- Ranger(kmer, concatenate)
  
  for (j in range) {
    #append new round to total
    kmer.colnames <- c(kmer.colnames, StatAppender(KmerGenerator(j), statistic))
  }
  
  kmer.colnames
}

################## calculate natural vectors ##################################
### input DNA sequence consisting of "a","c","g","t", updated on 11/08/2011
### output adjusted natural vector up to 40th dimension: nk, muk, D2k, E2k, E3k, ..., E8k
### output: vector -- 40-dim natural vector
###         length -- length of sequence
###         other  -- number of letters other than a,c,g,t
###         letter -- distinct letters in the sequence
CalculateVec <- function(kmer.seq, statistic = 3,
                         kmer = 3, method = "Sufficient") {
  
  seq.length <- length(kmer.seq)
  ans <- rep(0, 4 ^ kmer * statistic)  
  kmers <- KmerGenerator(k)
  
  #This is used to ensure there are only A, C, G, T.
  tbl <- table(c(kmer.seq, kmers)) - 1
  
  #There aren't/are other letters
  if (length(tbl) == length(kmers)) {
    if (sum(names(tbl) != kmers) == 0) {
      ans <- CalcRegLetters(kmer.seq, tbl, statistic, kmers, method)
    } else {
      ans <- CalcAmbigLetters(kmer.seq, statistic, kmers, method)
    }
  } else {
    if (sum(kmer.seq == "") > 0) {
      ans <- CalcEmptyLetters(tbl, statistic, kmers)
    } else {
      ans <- CalcAmbigLetters(kmer.seq, statistic, kmers, method)
    }
  }
  
  list(vector = ans, length = seq.length);
}

####################################################
ColumnFinder <- function(k, statistic, length.vec, concatenate) {
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
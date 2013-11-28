\name{Vectorizer}
\alias{Vectorizer}
\title{Calculate a generalized vector from a DNA sequence.}
\usage{
  Vectorizer(dna.seq, kmer = 3, statistic = 3,
    concatenate = TRUE, confirm.dna.seq = FALSE,
    verbose = 0)
}
\arguments{
  \item{dna.seq}{an upper or lower-case character vector
  consisting of A, C, T, G, or any of the IUPAC ambiguous
  nucleotide letters}

  \item{kmer}{an integer that defines the kmer length to be
  calculated}

  \item{concatenate}{logical indicating whether the
  statistics for the kmer should be returned or if the 1st
  through kth kmer statistics should be returned}

  \item{confirm.dna.seq}{logical indicating that the
  \code{dna.seq} should be confirmed as a proper dna
  sequence; this will add to the run-time}

  \item{verbose}{integer taking the values 0, 1, or 2 for
  varying levels of verbosity}
}
\description{
  This function takes a DNA sequence in character string
  format; e.g. "GATTACA" and calculates alignment free
  statistics.
}
\examples{
# K-mer with k = 2
Vectorizer(GetSeqGbk("data/abalone.gbk", upper = T),
                      kmer = 2, statistic = 1, concatenate = F)
#      length       nAA        nAC        nAG        nAT       nCA        nCC
# [1,]  34952 0.1059197 0.05338903 0.06194386 0.08572001 0.0779663 0.04094303
#            nCG        nCT        nGA        nGC        nGG        nGT        nTA
# [1,] 0.0219164 0.06283082 0.05267374 0.05130039 0.03705187 0.05052788 0.07044148
#           nTC        nTG        nTT
[1,] 0.05799548 0.07064176 0.09873823

# The Natural Vector
abalone <- Vectorizer(GetSeqGbk("data/abalone.gbk", upper = T),
                      kmer = 1, statistic = 3)
abalone.natural.vector <- abalone[1]*abalone[2:13]
abalone.natural.vector
# [1] 10730.000  7118.000  6695.000 10409.000 17319.117 17868.755 17451.659
# [8] 17386.478  2946.577  2794.014  2914.679  2954.686#'
}
\keyword{vectorizer}


\name{Vectorizer}
\alias{Vectorizer}
\title{Calculate a generalized vector from a DNA sequence.}
\usage{
Vectorizer(dna.seq, kmer = 3, statistic = 3, concatenate = TRUE,
  confirm.dna.seq = FALSE, verbose = 0)
}
\arguments{
  \item{dna.seq}{an upper or lower-case character vector
  consisting of A, C, T, G, or any of the IUPAC ambiguous
  nucleotide letters}

  \item{kmer}{an integer specifying the kmer length to be
  calculated}

  \item{statistic}{an integer from 1 to 5 specifying the
  statistics to be calculated; 1 corresponds to frequency,
  2 to frequency and mean, 3 includes variance, 4 skewness,
  and 5 kurtosis}

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
Vectorizer("GATTACA", kmer = 2, statistic = 1, concatenate = FALSE)
# length nAA       nAC nAG       nAT       nCA nCC nCG nCT       nGA nGC nGG nGT       nTA nTC nTG

# The Natural Vector
gattaca <- Vectorizer("GATTACA", kmer = 1, statistic = 3)
gattaca.natural.vector <- gattaca[1]*gattaca[2:13]
gattaca.natural.vector
}
\keyword{vectorizer}


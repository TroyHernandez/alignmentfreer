\name{CompTransform}
\alias{CompTransform}
\title{Calculate the composition vector from a generalized vector sequence.}
\usage{
  CompTransform(vec, kmer = 3, statistic = 3)
}
\arguments{
  \item{vec}{a generalized vector sequence as derived from
  \code{Vectorizer} with kmer greater than 1.}

  \item{kmer}{an integer that defines the kmer length of
  vec.}

  \item{statistic}{an integer that defines the statistics
  present in vec.}
}
\description{
  This function takes a generalized vector sequence, as
  obtained from \code{Vectorizer}, and calculates the
  composition vector.
}
\keyword{vectorizer}


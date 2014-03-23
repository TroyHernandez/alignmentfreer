\name{gbk}
\alias{gbk}
\title{Obtain the DNA sequence, and other information, from a gbk file.}
\usage{
gbk(path, upper = FALSE, phylo = NA)
}
\arguments{
  \item{path}{the path to a gbk file}

  \item{upper}{logical indicating if lower-case nucleotides
  should be converted to upper-case}

  \item{phylo}{string or \code{NA} indicating if
  phylogenetic information should be returned. Currently,
  only \code{"virus"} option is working.}
}
\description{
This function extracts a DNA sequence from a gbk file in
addition to other information; e.g. gi number, accession
number, number of base pairs, and organism name.  A
\code{gbk} class is output.
}
\keyword{vectorizer}


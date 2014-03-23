\name{fasta}
\alias{fasta}
\title{Obtain the DNA sequence, and other information, from a fasta file.}
\usage{
fasta(path, lower = FALSE)
}
\arguments{
  \item{path}{the path to a gbk file}

  \item{lower}{logical indicating if upper-case nucleotides
  should be converted to lower-case}
}
\description{
This function extracts a DNA sequence from a fasta file in
addition to other information; e.g. gi number, accession
number, number of base pairs, and organism name.  A
\code{fasta} class is output.
}
\keyword{vectorizer}


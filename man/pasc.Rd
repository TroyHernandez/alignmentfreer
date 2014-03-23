\name{pasc}
\alias{pasc}
\title{Convert a Pasc similarity file into a matrix.}
\usage{
pasc(path, similarity = TRUE)
}
\arguments{
  \item{path}{the path to a pasc similarity file}

  \item{similarity}{logical determining if similarity
  matrix or dissimilarity matrix should be returned}
}
\description{
This function converts a similarity file as obtained from
the PASC tool and converts it into a matrix and a
corresponding data frame.
}
\references{
http://www.ncbi.nlm.nih.gov/sutils/pasc/viridty.cgi
}
\keyword{vectorizer}


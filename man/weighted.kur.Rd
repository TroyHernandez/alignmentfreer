\name{weighted.kur}
\alias{weighted.kur}
\title{Compute a weighted kurtosis.}
\usage{
weighted.kur(x, w, na.rm = FALSE)
}
\arguments{
  \item{x}{an object containing the values whose weighted
  kurtosis is to be computed.}

  \item{w}{a numerical vector of weights the same length as
  x giving the weights to use for elements of x.}

  \item{na.rm}{a logical value indicating whether NA values
  in x should be stripped before the computation proceeds.}
}
\description{
Compute a weighted kurtosis.
}
\examples{
## GPA from Siegel 1994
wt <- c(5,  5,  4,  1)/15
x <- c(3.7,3.3,3.5,2.8)
xk <- weighted.kur(x, wt)
}
\keyword{vectorizer}


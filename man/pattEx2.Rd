\name{pattEx2}
\alias{pattEx2}
\docType{data}
\title{
Collection of two 2D patterns
}
\description{
Collection of two 2D patterns.
}
\usage{data(pattEx2)}
\format{
  A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{\code{key}}{A numeric vector indicating which points belong to which patterns.}
    \item{\code{x}}{A numeric vector corresponding to the first dimension of the points in the patterns.}
    \item{\code{y}}{A numeric vector corresponding to the second dimension of the points in the patterns.}
  }
}
\examples{
data(pattEx2)

#===> obtain the patterns <===#
x <- pattEx2[pattEx2[,1] == 1,2:3]
y <- pattEx2[pattEx2[,1] == 2,2:3]

#===> compute distance (via IMA) <===#
(d <- stDist(x, y, 0.5))
summary(d)
plot(d)
}
\keyword{datasets}

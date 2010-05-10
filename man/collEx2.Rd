\name{collEx2}
\alias{collEx2}
\docType{data}
\title{
Collection of two dimensional patterns
}
\description{
A collection of two dimensional patterns.
}
\usage{data(collEx2)}
\format{
  A data frame with 134 observations on the following 3 variables.
  \describe{
    \item{\code{key}}{A numeric vector corresponding to the pattern numbers.}
    \item{\code{x}}{A numeric vector corresponding to the first dimension of the points.}
    \item{\code{y}}{A numeric vector corresponding to the second dimension of the points.}
  }
}
\examples{
data(collEx2)

#===> construct the collection <===#
ppc   <- ppColl(collEx2[,2:3], collEx2[,1])
plot(ppc)

#===> compute the prototype via use margPT <===#
protoMP <- ppPrototype(ppc, 0.1)
points(protoMP, pch=20, cex=3.5, col='#FF000088')
# both algorithms produce very similar results
}
\keyword{datasets}

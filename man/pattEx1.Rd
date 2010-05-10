\name{pattEx1}
\alias{pattEx1}
\docType{data}
\title{
Collection of 2 single dimensional patterns
}
\description{
Collection of 2 single dimensional point patterns.
}
\usage{data(pattEx1)}
\format{
  A data frame with 16 observations on the following 2 variables.
  \describe{
    \item{\code{key}}{A numeric vector indicating which points belong to which patterns.}
    \item{\code{pts}}{A numeric vector corresponding to the points in the patterns.}
  }
}
\examples{
data(pattEx1)

#===> obtain the patterns <===#
x <- pattEx1[pattEx1[,1] == 1,2]
y <- pattEx1[pattEx1[,1] == 2,2]
par(mfrow=1:2)

#===> compute distance via VP97 <===#
(dVP <- stDist(x, y, 0.5, alg='VP97'))
summary(dVP)
plot(dVP)

#===> compute distance via IMA <===#
(dIMA <- stDist(x, y, 0.5, alg='IMA'))
summary(dIMA)
plot(dIMA)

#===> compute distance via MSU <===#
(dMSU <- stDist(x, y, 0.5, alg='MSU'))
summary(dMSU)
}
\keyword{datasets}

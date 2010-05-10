\name{collEx1}
\alias{collEx1}
\docType{data}
\title{
Collection of one dimensional patterns
}
\description{
A collection of one dimensional patterns.
}
\usage{data(collEx1)}
\format{
  A data frame with 332 observations on the following 2 variables.
  \describe{
    \item{\code{key}}{A numeric vector corresponding to the pattern numbers.}
    \item{\code{pts}}{A numeric vector corresponding to points.}
  }
}
\source{
Each pattern was simulated from an inhomogeneous Poisson process.
}
\examples{
data(collEx1)

#===> construct the collection <===#
ppc   <- ppColl(collEx1[,2], collEx1[,1])
plot(ppc, ylim=c(-3, 28))

#===> compute the prototype using VP97 <===#
protoVP <- ppPrototype(ppc, 0.5, alg='VP97')
points(protoVP, pch=20, cex=2, col='#FF000088')

#===>compute the prototype using margPT <===#
protoMP <- ppPrototype(ppc, 0.5)
points(protoMP, pch=20, cex=2, col='#00FF0088', at=-1)
# both algorithms produce very similar results
}
\keyword{datasets}

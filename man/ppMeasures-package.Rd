\name{ppMeasures-package}
\alias{ppMeasures-package}
\alias{ppMeasures}
\docType{package}
\title{
Point pattern distances and prototypes.
}
\description{
The package focuses on distances and prototypes for point patterns. There are three algorithms provided to compute spike-time distance, and one of these algorithms is generalized to compute variations of spike-time distance. Multiple algorithms are also provided to estimate prototypes of collections of point patterns.
}
\details{
\tabular{ll}{
Package: \tab ppMeasures\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-05-01\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
DM Diez, KE Tranbarger Freier, and FP Schoenberg

Maintainer: DM Diez <david.m.diez@gmail.com>
}
\references{
J Victor and K Purpura (1997). Metric-space analysis of spike trains:  theory, algorithms and application. Network: Computation in Neural Systems vol. 8, pp. 127-164, 1997.

KE Tranbarger and FP Schoenberg (2010). On the computation and application of point process prototypes.  Open Applied Informatics Journal 4, 1-9. 

DM Diez, FP Schoenberg, and CD Woody (as of April 2010, in review). Analysis of Neuronal Responses to Stimuli in Cats Using Point Process Prototypes.
}
\keyword{ package }
\keyword{ point pattern }
\keyword{ point process }
\keyword{ point pattern distance }
\keyword{ point pattern prototype }
\seealso{
\code{\link{stDist}}, \code{\link{ppColl}}, \code{\link{ppPrototype}}
}
\examples{
data(pattEx2)
x <- pattEx2[pattEx2[,1] == 1,c(2,3)]
y <- pattEx2[pattEx2[,1] == 2,c(2,3)]
(hold2 <- stDist(x, y, 2))
summary(hold2)
plot(hold2)

data(collEx2)
(ppc2 <- ppColl(collEx2[,2:3], collEx2[,1]))
summary(ppc2)
plot(ppc2, pch=0.5)
hold2 <- ppPrototype(ppc2, pm=0.05)
points(hold2, pch=20, cex=3, col='#FF000088')

#===> reproducing results from Diez et al. (in review) <===#
# results differ slightly due to default use of "margPT" algorithm
# and that cats are not weighted equally in the below analysis
par(mfrow=c(3,1))

data(neurNaive)
(ppc3 <- ppColl(neurNaive[,2], neurNaive[,1], nMissing=))
summary(ppc3)
plot(ppc3, cex=0.5)
hold3 <- ppPrototype(ppc3, pm=0.0106, bypassCheck=TRUE)
points(hold3, pch=20, cex=3, col='#FF000044')

data(neurCond)
(ppc4 <- ppColl(neurCond[,2], neurCond[,1]))
summary(ppc4)
plot(ppc4, cex=0.5)
hold4 <- ppPrototype(ppc4, pm=0.0075, bypassCheck=TRUE)
points(hold4, pch=20, cex=3, col='#FF000044')

data(neurBd)
(ppc5 <- ppColl(neurBd[,2], neurBd[,1]))
summary(ppc5)
plot(ppc5, cex=0.5)
hold5 <- ppPrototype(ppc5, pm=0.0078, bypassCheck=TRUE)
points(hold5, pch=20, cex=3, col='#FF000044')
}

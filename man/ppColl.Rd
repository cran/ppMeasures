\name{ppColl}
\alias{ppColl}
\alias{print.ppColl}
\alias{summary.ppColl}
\alias{plot.ppColl}
\alias{points.ppColl}
\title{
Creating a point pattern collection object
}
\description{
Prepare a point pattern collection for use in the prototype function, \code{\link{ppPrototype}}.
}
\usage{
ppColl(pts, key, wts = NULL, wtsKey = NULL, nMissing = 0)
}
\arguments{
  \item{pts}{A data matrix where rows represent points and columns dimensions.}
  \item{key}{A vector designating which points (rows of \code{pts}) belong to which patterns. Patterns without points should also have weights assigned.}
  \item{wts}{A numerical vector to specify weights for the patterns.}
  \item{wtsKey}{A vector. If \code{wts} are specified, then also specify the keys that those weights correspond to using \code{wtsKey}. This vector should have length equal to that of \code{wts} and include all unique values in \code{key}. Patterns without points should have a name assigned.}
  \item{nMissing}{The number of patterns with no points. This argument is ignored if weights are specified, in which case patterns without points should be assigned a weight in \code{wts} and a pattern name in \code{wtsKey}.}
}
\value{
An object of class \code{"ppColl"} is a list containing the following components:

\item{points }{A data matrix of the points.}
\item{pattID }{Original key proposed by the user.}
\item{key }{Key for identifying which points belong to which pattern.}
\item{keyMax }{The number of patterns.}
\item{wts }{Weights for each pattern.}
\item{dim }{Dimension of data.}
\item{maxObs }{The maximum number of observations in any of the patterns.}
\item{obs }{The number of observations in each pattern.}
}
\author{
DM Diez
}
\seealso{
\code{\link{ppPrototype}}, \code{\link{stDist}}
}
\examples{
data(collEx1)
(ppc1 <- ppColl(collEx1[,2], collEx1[,1]))
summary(ppc1)
plot(ppc1, ylim=c(-3, ppc1$keyMax+3))
hold1 <- ppPrototype(ppc1, pm=0.5)
points(hold1, pch=20, cex=3, col='#FF000088')

data(collEx2)
(ppc2 <- ppColl(collEx2[,2:3], collEx2[,1]))
summary(ppc2)
plot(ppc2, cex=0.5)
hold2 <- ppPrototype(ppc2, pm=0.05)
points(hold2, pch=20, cex=3, col='#FF000088')
}
\keyword{ Point pattern collection }
\keyword{ Point pattern prototype }

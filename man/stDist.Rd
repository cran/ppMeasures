\name{stDist}
\alias{stDist}
\alias{print.stDist}
\alias{summary.stDist}
\alias{plot.stDist}
\title{
Spike-time distance computation
}
\description{
This function computes spike-time distance and its variants.
}
\usage{
stDist(x, y, pm, pa = 1, pd = 1, alg = c("default", "IMA", "VP97", "MSU"), euclid = NULL,
       lossOrder = 1, maxBranch = 4, eps = 10^-10, bypassCheck = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix representing a point pattern where rows represent points and columns dimensions. If there is only one dimension, \code{x} may be a vector.}
  \item{y}{A matrix representing a point pattern where rows represent points and columns dimensions. If there is only one dimension, \code{y} may be a vector.}
  \item{pm}{A number or numerical vector representing the moving penalty. If there are multiple dimensions, different moving penalties may be specified for each dimension by using a vector for \code{pm}.}
  \item{pa}{A number representing the addition penalty.}
  \item{pd}{A number representing the deletion penalty.}
  \item{alg}{Character string specifying which algorithm should be applied. \code{"IMA"} corresponds to an incremental matching algorithm for spike-time distance and its variants in any finite dimensional space. \code{"VP97"} corresponds to a slightly generalized algorithm from Victor and Purpura (1997) for spike-time distance. \code{"MSU"} corresponds to the Modified Single Unit algorithm for spike-time distance (Tranbarger and Schoenberg, 2010). By default, \code{"VP97"} will be used for one-dimensional data and \code{"IMA"} will be used for multiple dimensions. If variants of spike-time distance are used (\code{euclid} or \code{lossOrder} are adjusted), then the \code{"IMA"} technique will be applied regardless of the specified algorithm.}
  \item{euclid}{A vector of \code{1}s and \code{0}s. The vector should have one entry per dimension, and dimensions corresponding to a \code{1} are incorporated as a Euclidean dimension for Euclidean distance between points.}
  \item{lossOrder}{A number. The moving penalty is computed using the moving penalty and distance a point is moved. The \code{lossOrder} specifies the power of the distance moved in the penalty computation, e.g. for one non-Euclidean dimension, the cost of moving a point a distance \code{d} is \code{pm*d^lossOrder}.}
  \item{maxBranch}{A non-negative integer used as a parameter for accuracy when \code{alg="IMA"}. Values larger than 2 typically result in an exact computation of distance while using a \code{maxBranch} of \code{0} or \code{1} will reduce accuracy but results in a faster computation.}
  \item{eps}{If reductions in distance are not larger than this amount, the algorithm will stop.}
  \item{bypassCheck}{If patterns may be too large for the specified algorithm, an error will occur. To bypass this error and proceed anyways, specify \code{bypassCheck=TRUE}.}
}
\details{
This function computes the distance from \code{x} to \code{y}.
}
\value{
An object of class \code{"stDist"} is a list containing the following components:

\item{distance }{The computed distance.}
\item{x }{The first pattern.}
\item{y }{The second pattern.}
\item{pa }{Addition penalty.}
\item{pd }{Deletion penalty.}
\item{pm }{Moving penalty.}
\item{xtoy }{Under the optimal matching, \code{x[i,]} matches to \code{y[xtoy[i],]}. Not available when \code{alg="MSU"}.}
\item{ytox }{Under the optimal matching, \code{y[i,]} matches to \code{x[ytox[i],]}. Not available when \code{alg="MSU"}.}
\item{dim }{Number of dimensions.}
\item{eps }{The specified \code{eps}.}
\item{maxBranch }{The specified \code{maxBranch}.}
\item{lossOrder }{The specified \code{lossOrder}.}
\item{algorithm }{The algorithm utilized.}
\item{euclid }{The specified \code{euclid}.}
}

\references{
J Victor and K Purpura (1997). Metric-space analysis of spike trains:  theory, algorithms and application. Network: Computation in Neural Systems vol. 8, pp. 127-164, 1997.

KE Tranbarger and FP Schoenberg (2010). On the computation and application of point process prototypes.  Open Applied Informatics Journal 4, 1-9. 

DM Diez, FP Schoenberg, and CD Woody (as of May 2010, in review). Analysis of Neuronal Responses to Stimuli in Cats Using Point Process Prototypes.
}

\author{
DM Diez, KE Tranbarger Freier, and FP Schoenberg
}

\seealso{
\code{\link{ppColl}}, \code{\link{ppPrototype}}
}
\examples{
data(pattEx1)
x <- pattEx1[pattEx1[,1] == 1,2]
y <- pattEx1[pattEx1[,1] == 2,2]
(hold1  <- stDist(x, y, 4))
(hold11 <- stDist(x, y, 4, alg='MSU'))
summary(hold1)
summary(hold11)
plot(hold1)

data(pattEx2)
x <- pattEx2[pattEx1[,1] == 1,c(2,3)]
y <- pattEx2[pattEx1[,1] == 2,c(2,3)]
(hold2 <- stDist(x, y, 2))
summary(hold2)
plot(hold2)
}
\keyword{ spike-time distance }
\keyword{ point pattern }

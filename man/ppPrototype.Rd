\name{ppPrototype}
\alias{ppPrototype}
\alias{print.ppPrototype}
\alias{summary.ppPrototype}
\alias{plot.ppPrototype}
\alias{points.ppPrototype}
\title{
Point pattern prototype estimation
}
\description{
This function estimates the prototype of a collection of point patterns.
}
\usage{
ppPrototype(ppcoll, pm, pa = 1, pd = 1, alg = c("margPT", "kernPT", "VP97", "MSU"),
            costEps = 10^(-6), posEps = 10^(-6), ppd = 4, space = 0.1, euclid = NULL,
            lossOrder = 1, maxBranch = 4, bypassCheck = FALSE)
}
\arguments{
  \item{ppcoll}{A collection of point patterns, output from the function \code{\link{ppColl}}.}
  \item{pm}{A number or numerical vector representing the moving penalty. If there are multiple dimensions, different moving penalties may be specified for each dimension by using a vector for \code{pm}.}
  \item{pa}{A number representing the addition penalty.}
  \item{pd}{A number representing the deletion penalty.}
  \item{alg}{Character string specifying the algorithm. \code{"margPT"} corresponds to a marginal prototype technique tied that utilizes the Incremental Matching Algorithm, which is suitable for multiple dimensions. Other (one-dimensional) methods check all possible prototype point locations and utilize different algorithms for the distance: \code{"kernPT"} for kernel smoothing, \code{"VP97"} for the algorithm presented in Victor and Purpura (1997), or \code{"MSU"} for the Modified Single Unit Algorithm.}
  \item{costEps}{If this reduction in cost is not achieved, the algorithm will stop.}
  \item{posEps}{A number or numerical vector for precision control for identification of prototype points. For instance, if the algorithm converges to a region of diameter smaller than \code{posEps}, it will complete its search for the prototype point. If working in multidimensional space and a numerical vector is used, the value in position i of \code{posEps} corresponds to the precision control for dimension i.}
  \item{ppd}{An integer or numerical vector of integers utilized in \code{alg="margPT"} to specify the number of points per dimension to examine at a single time. Each value must be at least \code{2}. If a vector, element i corresponds to the number of points to examine in dimension i at a single time.}
  \item{space}{A value greater than 0 but less than 0.5, which governs the spread of the points in each dimension when using \code{alg="margPT"}. Just as with \code{ppd}, space can be input as a vector where the i-th value corresponds to dimension i. Restriction on \code{space} and \code{ppd}: for any dimension, their product must be less than 0.5.}
  \item{euclid}{A vector of \code{1}s and \code{0}s. The vector should have one entry per dimension, and dimensions corresponding to a \code{1} are incorporated as a Euclidean dimension for Euclidean distance between points. This argument is only utilized when \code{alg="margPT"}.}
  \item{lossOrder}{A number. The moving penalty is computed using the moving penalty and distance a point is moved. The \code{lossOrder} specifies the power of the distance moved in the penalty computation, e.g. for one non-Euclidean dimension, the cost of moving a point a distance \code{d} is \code{pm*d^lossOrder}. This argument is only utilized when \code{alg="margPT"}.}
  \item{maxBranch}{A non-negative integer used as a parameter for accuracy when \code{alg="margPT"}. Values larger than 2 typically result in an exact computation of distance while using a \code{maxBranch} of \code{0} or \code{1} will reduce accuracy but results in a faster computation.}
  \item{bypassCheck}{If patterns may be too large for the specified algorithm, an error will occur. To bypass this error and proceed anyways, specify \code{bypassCheck=TRUE}.}
}
\details{
Doing a test with a smaller data set is recommended since large data sets may require a significant amount of time for computation. Computation time is dependent not only on the number of patterns or the number of points in the patterns, but also dimension, \code{ppd}, and other settings and collection characteristics.
}
\value{
A list containing the following items:

 \item{algorithm }{The algorithm utilized.}
 \item{ppColl }{The collection for which the prototype was found.}
 \item{ppColl }{The collection for which the prototype was found.}
 \item{prototype }{Matrix representing the prototype.}
 \item{match }{How points in the collection are matched to the prototype.}
 \item{cost }{Total cost: weighted sum of the distances from the prototype to each pattern.}
 \item{pa }{Addition penalty.}
 \item{pa }{Deletion penalty.}
 \item{pa }{Moving penalty.}
 \item{ppd }{Points per dimension, if utilized the \code{"margPT"} algorithm.}
 \item{space }{Spacing parameter of prototype search, if utilized the \code{"margPT"} algorithm.}
 \item{dim }{Dimension of the space.}
 \item{costEps }{\code{costEps} parameter.}
 \item{posEps }{\code{posEps} parameter.}
 \item{maxBranch }{\code{maxBranch} parameter.}
 \item{lossOrder }{\code{lossOrder} parameter.}
 \item{euclid }{\code{euclid} parameter.}
 \item{dist }{The distances from the prototype to each pattern in the collection.}
}
\references{
KE Tranbarger and FP Schoenberg (2010). On the computation and application of point process prototypes.  Open Applied Informatics Journal 4, 1-9. 

Diez DM, Schoenberg FP, and Woody CD (2012). Algorithms for computing spike time
distance and point process prototypes with application to feline neuronal responses to
acoustic stimuli. Journal of Neuroscience Methods 203(1):186-192.
}
\author{
DM Diez
}

\seealso{
\code{\link{ppColl}}, \code{\link{stDist}}
}
\examples{
#===> one dimensional example <===#
data(collEx1)
ppc   <- ppColl(collEx1[,2], collEx1[,1])
plot(ppc, ylim=c(-3, 28))
protoVP <- ppPrototype(ppc, 0.5, alg='VP97')
summary(protoVP)
points(protoVP, pch=20, cex=2, col='#FF000088')
protoMP <- ppPrototype(ppc, 0.5)
points(protoMP, pch=20, cex=2, col='#00FF0088', at=-1)
protoKP <- ppPrototype(ppc, 0.5, alg='kernPT', ppd=500)
points(protoKP, pch=20, cex=2, col='#0000FF88', at=-2)
protoSP <- ppPrototype(ppc, 0.5, alg='MSU', ppd=500, bypassCheck=TRUE)
points(protoSP, pch=20, cex=2, col='#FFFF0088', at=-3)

#===> two dimensional example <===#
data(collEx2)
ppc   <- ppColl(collEx2[,2:3], collEx2[,1])
plot(ppc)
protoMP <- ppPrototype(ppc, 0.1)
points(protoMP, pch=20, cex=3.5, col='#FF000088')
}
\keyword{ point pattern prototypes }
\keyword{ point pattern collection }
\keyword{ spike-time distance }

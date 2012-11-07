\name{neurCond}
\alias{neurCond}
\docType{data}
\title{
Collection of cell prototypes from conditioned cats
}
\description{
Collection of cell prototypes from conditioned cats. This data originates from an experiment testing the effect of conditioning and backwards conditioning on the neural responses of cats (Woody et al., 2000). This data pertains to the cats following stimulus conditioning and their neuron cells' prototypical responses to a click stimulus at time 0 ms that were observed between -320 ms and 1280 ms. Originally each cell had multiple runs, and these have been simplified to cell prototypes in this data set.
}
\usage{data(neurCond)}
\format{
  A data frame with 2855 observations on the following 3 variables.
  \describe{
    \item{\code{key}}{A numeric vector corresponding to the cell numbers.}
    \item{\code{pts}}{A numeric vector corresponding to neuronal responses.}
    \item{\code{cat}}{A numeric vector corresponding to the cat.}
  }
}
\details{
Each cell originally had several observed patterns. To simplify this complex collection, a prototype was found for each cell.

An analysis incorporating the weights of the cats is included in the examples (Diez et al., 2010). Thirty-five cells had zero points in their prototype (6 for cat \code{1}, 4 for \code{2}, 15 for \code{3}, 5 for \code{4}, 5 for \code{5}).
}
\source{
CD Woody, E Zotova, and E Gruen (2000). Multiple representations of information in the primary auditory cortex of cats: I. stability and change in slow components of unit activity after conditioning with a click conditioned stimulus. Brain Research, 868, 56-65.
}
\references{
Diez DM, Schoenberg FP, and Woody CD (2012). Algorithms for computing spike time
distance and point process prototypes with application to feline neuronal responses to
acoustic stimuli. Journal of Neuroscience Methods 203(1):186-192.
}
\examples{
data(neurCond)

#===> construct weights <===#
# wts found separately
counts <- c(6,4,15,5,5)
wtsKey <- 1:256
# cat information from original data
cats   <- c(rep(1,29), rep(2,55), rep(3,31), rep(4,16),
		rep(3,54), rep(5,44), rep(3,27))
for(i in 1:5){
	these <- unique(neurCond[neurCond[,3] == i,1])
	counts[i] <- counts[i] + length(these)
}
sum(counts)
wts <- 1/counts[cats]

#===> setup collection <===#
pts <- neurCond[,2]
key <- neurCond[,1]
ppc <- ppColl(pts, key, wts, wtsKey)
# did not need to list nMissing since implied through wts/wtsKey

#===> identify prototype <===#
(protoKP <- ppPrototype(ppc, 0.0075, pa=1.1125, pd=0.8875, alg='kernPT', ppd=5000, bypassCheck=TRUE))
plot(ppc)
points(protoKP, col=4, cex=1.5, at=-10)
}
\keyword{datasets}

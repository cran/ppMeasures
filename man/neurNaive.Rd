\name{neurNaive}
\alias{neurNaive}
\docType{data}
\title{
Collection of cell prototypes from naive cats
}
\description{
Collection of cell prototypes from naive cats. This data originates from an experiment testing the effect of conditioning and backwards conditioning on the neural responses of cats (Woody et al., 2000). This data pertains to naive cats (prior to stimulus conditioning) and their neuron cells' prototypical responses to a click at time 0 ms that were observed between -320 ms and 1280 ms. Originally each cell had multiple runs, and these have been simplified to cell prototypes in this data set.
}
\usage{data(neurNaive)}
\format{
  A data frame with 3177 observations on the following 3 variables.
  \describe{
    \item{\code{key}}{A numeric vector corresponding to the cell numbers.}
    \item{\code{pts}}{A numeric vector corresponding to neuronal responses.}
    \item{\code{cat}}{A numeric vector corresponding to the cat.}
  }
}
\details{
Each cell originally had several observed patterns. To simplify this complex collection, a prototype was found for each cell.

An analysis incorporating the weights of the cats is included in the examples (Diez et al., 2010). Twelve cells had zero points in their prototype (3 for cat \code{1}, 1 for \code{2}, 3 for \code{3}, 5 for \code{4}).
}
\source{
CD Woody, E Zotova, and E Gruen (2000). Multiple representations of information in the primary auditory cortex of cats: I. stability and change in slow components of unit activity after conditioning with a click conditioned stimulus. Brain Research, 868, 56-65.
}
\references{
DM Diez, FP Schoenberg, and CD Woody (as of April 2010, in review). Analysis of Neuronal Responses to Stimuli in Cats Using Point Process Prototypes.
}
\examples{
data(neurNaive)

#===> construct weights <===#
# wts found separately
counts <- c(3,1,3,5)
wtsKey <- 1:268
# cat information from original data
cats   <- c(rep(1,66), rep(2,39), rep(3,80), rep(4,83))
for(i in 1:4){
	these <- unique(neurNaive[neurNaive[,3] == i,1])
	counts[i] <- counts[i] + length(these)
}
sum(counts)
wts <- 1/counts[cats]

#===> setup collection <===#
pts <- neurNaive[,2]
key <- neurNaive[,1]
ppc <- ppColl(pts, key, wts, wtsKey)

#===> identify prototype <===#
(protoKP <- ppPrototype(ppc, 0.0106, pa=1.1, pd=0.9, alg='kernPT', ppd=5000, bypassCheck=TRUE))
plot(ppc)
points(protoKP, col=4, cex=1.5, at=-10)
}
\keyword{datasets}

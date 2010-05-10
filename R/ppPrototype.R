ppPrototype <-
function(ppcoll, pm, pa=1, pd=1, alg=c('margPT', 'kernPT', 'VP97', 'MSU'), costEps=10^(-6), posEps=10^(-6), ppd=4, space=0.1, euclid=NULL, lossOrder=1, maxBranch=4, bypassCheck=FALSE){
	tR <- list()
	class(tR)    <- 'ppPrototype'
	tR$algorithm <- NA
	tR$ppColl    <- NA
	tR$prototype <- NA
	tR$match <- NA
	tR$cost  <- NA
	tR$pa    <- pa
	tR$pd    <- pd
	tR$pm    <- pm
	tR$ppd   <- ppd
	tR$space <- space
	tR$dim   <- NA
	tR$dist  <- c()
	tR$costEps   <- costEps
	tR$posEps    <- NA
	tR$maxBranch <- maxBranch
	tR$euclid    <- NA
	tR$lossOrder <- lossOrder

#	suppress for now, had argument kernProto=FALSE previously
#	tR$kernProto <- sum(kernProto)
#	if(tR$kernProto == 1 & tR$maxBranch != 0){
#		cat('Warning: Because maxBranch was not 0, kernProto was reset to 0.\n')
#	}
	
	if(class(ppcoll) != 'ppColl'){
		stop('Argument "ppcoll" must be of class "ppColl".\n')
	} else {
		tR$ppColl <- ppcoll
		tR$dim    <- ppcoll$dim
	}
	npts <- dim(ppcoll$points)[1]
	maxObs <- sum(ppcoll$key==1)
	for(i in 2:ppcoll$keyMax){
		temp <- sum(ppcoll$key==i)
		if(temp > maxObs){
			maxObs <- temp
		}
	}
	if(length(tR$pm) != tR$dim){
		tR$pm <- rep(tR$pm, tR$dim)
	}
	if(length(tR$ppd) < tR$dim){
		tR$ppd <- rep(tR$ppd, tR$dim)
	}
	if(length(tR$space) < tR$dim){
		tR$space <- rep(tR$space, tR$dim)
	}
	PT  <- rep(-1, maxObs*tR$dim)
	nPT <- 0
	PTcost   <- -1
	tR$match <- rep(-1, npts)
	if(is.null(euclid)[1]){
		tR$euclid <- c()
	} else {
		if(is.numeric(euclid)){
			tR$euclid <- euclid
		} else {
			if(length(euclid) == tR$dim){
				tR$euclid <- which(euclid)
			} else {
				stop('Argument "euclid" is not recognized.\n')
			}
		}
	}
	if(length(posEps) < tR$dim){
		tR$posEps <- rep(posEps, tR$dim)
	} else {
		tR$posEps <- posEps
	}
	temp <- rep(0, tR$dim)
	temp[tR$euclid] <- 1
	tR$euclid <- temp
	if(alg[1] == 'kernPT'){
		if(tR$dim == 1){
			kernPT <- TRUE
			margPT <- FALSE
		} else {
			stop('The kernel smoothing approach is only permitted in a single dimension in this current release.\n')
		}
	} else if(alg[1] == 'VP97' & tR$dim == 1){
		VP97   <- TRUE
		margPT <- FALSE
		kernPT <- FALSE
	} else if(alg[1] == 'MSU' & tR$dim == 1){
		MSU    <- TRUE
		margPT <- FALSE
		kernPT <- FALSE
		VP97   <- FALSE
	} else {
		margPT <- TRUE
	}
	if(margPT){
		temp <- rep(100, floor(sqrt(tR$ppColl$keyMax)))
		temp <- c(temp, tR$ppColl$obs)
		c0 <- quantile(temp, pa/(pa+pd)) > 25
		c1 <- tR$ppColl$maxObs > 100
		c2 <- sum(tR$ppColl$obs > 80) > 5
		c2 <- c2 && (sum(tR$ppColl$obs > 60) > 10)
		c3 <- sum(tR$ppColl$obs > 35) > 15
		c4 <- sum(tR$ppColl$obs > 15) > 60
		c5 <- sum(tR$ppColl$obs > 10) > 120
		c6 <- sum(tR$ppColl$obs > 8) > 250
		cc <- (c0 || c1 || c2 || c3 || c4 || c5 || c6)
		if(!bypassCheck && cc){
			stop('Running "ppPrototype" on this collection under the "margPT" algorithm may freeze or crash your R session. If you wish to proceed anyways, use the "bypassCheck" argument.\n')
		}
		if(any(ppd*space >= 0.5)){
			stop('The product of the arguments "ppd" and "space" must be less than 0.5 (for each dimension).\n')
		}
		tR$algorithm <- 'margPT'
		cOut  <- .C("margPT",
			as.double(ppcoll$points),  # 1
			as.integer(npts),          # 2
			as.integer(ppcoll$key),    # 3
			as.integer(tR$dim),        # 4
			as.integer(ppcoll$keyMax), # 5
			as.integer(maxObs), # 6
			as.double(tR$pm),   # 7
			as.double(tR$pa),   # 8
			as.double(tR$pd),   # 9
			as.double(ppcoll$wts),    # 10
			as.integer(tR$maxBranch), # 11
			as.integer(tR$ppd),       # 12
			as.double(tR$space),      # 13
			as.double(PT),     # 14
			as.integer(nPT),   # 15
			as.double(PTcost), # 16
			as.double(tR$costEps), # 17
			as.double(tR$posEps),  # 18
			as.integer(tR$match),  # 19
			as.integer(tR$euclid), # 20
			as.double(tR$lossOrder),    # 21
			as.integer(tR$useKernOnly), # 22 ignored
			as.integer(0), # 23 kernProto uses other fcn
			PACKAGE="ppMeasures")
		if(cOut[[15]] > 0){
			PT <- matrix(cOut[[14]][1:(tR$dim*cOut[[15]])],
				ncol=tR$dim, byrow=TRUE)
		} else {
			PT <- NULL
		}
		tR$cost  <- cOut[[16]]
		tR$match <- cOut[[19]]+1
	} else if(kernPT){
		# need PT, tR$cost, tR$match
		tR$algorithm <- 'kernPT'
		c1 <- tR$ppColl$maxObs > 120
		c2 <- length(tR$ppColl$key)/tR$ppColl$keyMax > 40
		c3 <- length(tR$ppColl$key)^tR$dim > 5000
		c4 <- sum(tR$ppColl$obs > 25) > 100
		cc <- (c1 || c2 || c3 || c4)
		if(!bypassCheck && cc){
			stop('Running "ppPrototype" on this collection under the "kernPT" algorithm may freeze or crash your R session. If you wish to proceed anyways, use the "bypassCheck" argument.\n')
		}
		cOut <- .C("kernProto",
			as.double(ppcoll$points), # 1
			as.integer(npts),         # 2
			as.integer(ppcoll$key),   # 3
			as.integer(tR$match),     # 4
			as.integer(ppcoll$keyMax), # 5
			as.integer(tR$dim),        # 6
			as.double(tR$pm), # 7
			as.double(tR$pa), # 8
			as.double(tR$pd), # 9
			as.double(ppcoll$wts), # 10
			as.integer(tR$euclid),   # 11
			as.double(tR$lossOrder), # 12
			as.double(PT),   # 13
			as.integer(nPT), # 14
			as.double(tR$posEps), # 15
			as.integer(tR$ppd),   # 16
			as.integer(prod(tR$ppd)), # 17
			as.double(0),             # 18 cost
			PACKAGE="ppMeasures")
		if(cOut[[14]] > 0){
			PT <- matrix(cOut[[13]][1:(tR$dim*cOut[[14]])],
				ncol=tR$dim, byrow=TRUE)
		} else {
			PT <- NULL
		}
		tR$cost  <- cOut[[18]]
		tR$match <- cOut[[4]]+1
	} else if(VP97){
		tR$algorithm <- 'VP97'
		temp1 <- length(tR$ppColl$key)
		temp2 <- quantile(tR$ppColl$obs, pa/(pa+pd))
		c1 <- (temp1*temp2)^2 > 2.5*10^8
		if(!bypassCheck && c1){
			stop('Running "ppPrototype" on this collection under the "VP97" algorithm may freeze or crash your R session. If you wish to proceed anyways, use the "bypassCheck" argument.\n')
		}
		cOut  <- .C("vpProto",
			as.double(ppcoll$points),  # 1
			as.integer(npts),          # 2
			as.integer(ppcoll$key),    # 3
			as.integer(ppcoll$keyMax), # 4
			as.integer(maxObs), # 5
			as.double(tR$pm),   # 6
			as.double(tR$pa),   # 7
			as.double(tR$pd),   # 8
			as.double(ppcoll$wts),    # 9
			as.double(PT),            # 10
			as.integer(nPT),          # 11
			as.double(PTcost),        # 12
			as.double(tR$costEps), # 13
			as.double(tR$posEps),  # 14
			as.integer(tR$match),    # 15
			as.double(tR$lossOrder), # 16 ignored
			PACKAGE="ppMeasures")
		if(cOut[[11]] > 0){
			PT <- matrix(cOut[[10]][1:(tR$dim*cOut[[11]])],
				ncol=tR$dim, byrow=TRUE)
		} else {
			PT <- vector()
		}
		tR$cost  <- cOut[[12]]
		tR$match <- cOut[[15]] + 1 # ZZQ make sure +1 is right!
	} else { # MSU
		tR$algorithm <- 'MSU'
		c1 <- tR$ppColl$maxObs > 23
		c2 <- length(tR$ppColl$key)/tR$ppColl$keyMax > 10
		c3 <- length(tR$ppColl$key)^tR$dim > 200
		if(!bypassCheck && (c1 || c2 || c3)){
			stop('Running "ppPrototype" on this collection under the "MSU" algorithm may freeze or crash your R session. If you wish to proceed anyways, use the "bypassCheck" argument.\n')
		}
		# need PT, tR$cost, tR$match
		cOut  <- .C("msuProto",
			as.double(ppcoll$points),  # 1
			as.integer(npts),          # 2
			as.integer(ppcoll$key),    # 3
			as.integer(ppcoll$keyMax), # 4
			as.integer(maxObs), # 5
			as.double(tR$pm),   # 6
			as.double(tR$pa),   # 7
			as.double(tR$pd),   # 8
			as.double(ppcoll$wts),    # 9
			as.double(PT),            # 10
			as.integer(nPT),          # 11
			as.double(PTcost),        # 12
			as.double(tR$costEps), # 13
			as.double(tR$posEps),  # 14
			as.integer(tR$match),    # 15
			as.double(tR$lossOrder), # 16 ignored
			PACKAGE="ppMeasures")
		if(cOut[[11]] > 0){
			PT <- matrix(cOut[[10]][1:(tR$dim*cOut[[11]])],
				ncol=tR$dim, byrow=TRUE)
		} else {
			PT <- vector()
		}
		tR$cost  <- cOut[[12]]
	}
	
	Alg <- tR$algorithm
	if(Alg %in% c('margPT', 'kernPT')){
		Alg <- 'IMA'
	}
	
	tR$prototype <- PT
	tR$match <- 0
	if(length(tR$prototype) > 0){
		for(i in 1:ppcoll$keyMax){
			these <- ppcoll$key == i
			if(sum(these) > 0 && length(tR$prototype) > 0){
				temp <- stDist(x=tR$prototype,
					y=ppcoll$points[these,],
					pm=tR$pm, pa=tR$pa, pd=tR$pd,
					eps=tR$costEps,
					euclid=tR$euclid,
					lossOrder=tR$lossOrder,
					maxBranch=tR$maxBranch,
					alg=Alg,
					bypassCheck=bypassCheck)
				tR$dist[i] <- temp$distance
				if(Alg != 'MSU'){
					tR$match[these] <- temp$ytox
				} else {
					tR$match[these] <- NA
				}
			} else {
				tR$dist[i] <- pd*dim(tR$prototype)[1]
				if(sum(these) > 0){
					tR$match[these] <- 0
				}
			}
		}
	} else {
		for(i in 1:ppcoll$keyMax){
			temp <- (ppcoll$key == i)
			tR$dist[i] <- pa*sum(temp)
		}
	}
	return(tR)
}


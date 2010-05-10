stDist <-
function(x, y, pm, pa=1, pd=1, alg=c('default', 'IMA', 'VP97', 'MSU'), euclid=NULL, lossOrder=1, maxBranch=4, eps=10^-10, bypassCheck=FALSE){
	tR <- list()
	class(tR)   <- 'stDist'
	tR$distance <- NA
	tR$x  <- NA
	tR$y  <- NA
	tR$pa <- pa
	tR$pd <- pd
	tR$pm <- NA
	tR$xtoy <- NA
	tR$ytox <- NA
	tR$dim  <- NA
	tR$eps  <- eps
	tR$maxBranch <- maxBranch
	tR$euclid    <- NA
	tR$lossOrder <- NA
	tR$algorithm <- alg[1]
	
	#if(length(y) > 0){
	#	if(is.na(y)[1]){
	#		if(class(x) == 'stDist'){
	#			stdist <- x
	#			x <- stdist$x
	#			y <- stdist$y
	#		} else {
	#			stop('Arguments "y" is required unless "x" is of class "stDist".\n')
	#		}
	#	}
	#}
	
	if(is.data.frame(x)){
		x <- as.matrix(x)
	}
	if(is.vector(x)){
		tR$x <- matrix(x, ncol=1)
	} else if(length(x) == 0){
		tR$x <- vector()
		nx   <- 0
	} else {
		tR$x <- x
	}
	if(is.data.frame(y)){
		y <- as.matrix(y)
	}
	if(is.vector(y)){
		tR$y <- matrix(y, ncol=1)
		ny   <- 0
	} else if(length(y) == 0){
		tR$y <- vector()
		ny   <- 0
	} else {
		tR$y <- y
	}

	if(length(x) == 0 & length(y) == 0){
		tR$dim = 1
	} else if(length(y) > 0){
		tR$dim = dim(tR$y)[2]
	} else {
		tR$dim = dim(tR$x)[2]
	}

	if(length(x) > 0){
		nx   <- dim(tR$x)[1]
	}
	if(length(y) > 0){
		ny   <- dim(tR$y)[1]
	}
	
	if(length(pm) == 1 & tR$dim > 1){
		tR$pm <- rep(pm, tR$dim)
	} else {
		tR$pm <- pm
	}
	tR$xtoy <- rep(-1, nx)
	tR$ytox <- rep(-1, ny)
	cost <- tR$pd*nx + tR$pa*ny
	
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
	temp <- rep(0, tR$dim)
	temp[tR$euclid] <- 1
	tR$euclid <- temp
	
	if(!is.numeric(lossOrder)){
		stop('Argument "lossOrder" must be a number greater than zero.\n')
	}
	if(!(lossOrder > 0)){
		stop('Argument "lossOrder" must be a number greater than zero.\n')
	}
	tR$lossOrder <- lossOrder
	
	temp <- min(c(nx, ny))
	nB <- 2*max((nx-1:temp)*(ny-1:temp)*(1:temp))
	mustUseIMA <- !is.null(euclid[1]) || any(lossOrder != 1)
	MSU <- (tR$dim == 1 &
			tR$algorithm == 'MSU' &
			tR$pa == tR$pd & !mustUseIMA)
	VP97 <- (tR$dim == 1 &
			(tR$algorithm == 'default' |
				tR$algorithm == 'VP97') & !mustUseIMA)
	cc <- c(length(x), length(y))
	if(bypassCheck){
		# let it run!
	} else if(MSU){
		c1 <- any(cc > 29)
		c2 <- prod(cc) > 300
		if(c1 || c2){
			stop('Due to the size of the patterns, the "MSU" algorithm may freeze your R session. To suppress this warning and proceed anyways, either switch algorithms or utilize the "bypassCheck" argument.\n')
		}
	} else if(VP97){
		c1 <- prod(cc) > 500^2
		if(c1){
			stop('Due to the size of the patterns being used, using "stDist" may freeze or crash your R session. To suppress this warning and proceed anyways, you may utilize the "bypassCheck" argument.\n')
		}
	} else {
		c1 <- prod(cc) > 100^2
		c2 <- any(cc > 300)
		if(c1 || c2){
			if(tR$dim > 1){
				stop('Due to the size of the patterns being used, using "stDist" may freeze or crash your R session. To suppress this warning and proceed anyways, you may utilize the "bypassCheck" argument.\n')
			} else {
				stop('Due to the size of the patterns being used, using the "IMA" algorithm may freeze or crash your R session. You might instead try the "VP97" algorithm. To suppress this warning and proceed anyways, you may utilize the "bypassCheck" argument.\n')
			}
		}
	}
	if(MSU){
		tR$x <- matrix(sort(tR$x), ncol=1)
		tR$y <- matrix(sort(tR$y), ncol=1)
		cOut <- .C("msuDist",
			as.integer(nx),   # 1
			as.integer(ny),   # 2
			as.double(tR$x),  # 3
			as.double(tR$y),  # 4
			as.double(tR$pa), # 5
			as.double(tR$pm), # 6
			#as.double(tR$pd),
			as.double(cost),  # 7
			PACKAGE="ppMeasures")
		tR$distance <- cOut[[7]]
		tR$xtoy <- NA
		tR$ytox <- NA
#		cat('MSU\n')
	} else if(VP97){
		tR$x <- matrix(sort(tR$x), ncol=1)
		tR$y <- matrix(sort(tR$y), ncol=1)
		cOut <- .C("vpAlg",
			as.double(tR$x),  # 1
			as.double(tR$y),  # 2
			as.integer(nx),   # 3
			as.integer(ny),   # 4
			as.double(tR$pm), # 5
			as.double(tR$pa), # 6
			as.double(tR$pd), # 7
			as.integer(tR$xtoy), # 8
			as.integer(tR$ytox), # 9
			as.double(cost),     # 10
			PACKAGE="ppMeasures")
			
		tR$distance <- cOut[[10]]
		tR$xtoy <- cOut[[8]] + 1
		tR$ytox <- cOut[[9]] + 1
		tR$algorithm <- 'VP97'
#		cat('VP97\n')
	} else {
		cOut <- .C("stdist",
			as.double(tR$x), # 1
			as.double(tR$y), # 2
			as.integer(nx),  # 3
			as.integer(ny),  # 4
			as.integer(tR$dim), # 5
			as.double(tR$pm),   # 6
			as.double(tR$pa),   # 7
			as.double(tR$pd),   # 8
			as.double(eps),     # 9
			as.integer(tR$xtoy),      # 10
			as.integer(tR$ytox),      # 11
			as.integer(tR$maxBranch), # 12
			as.double(cost),          # 13
			as.integer(nB),           # 14
			as.integer(tR$euclid),    # 15
			as.double(tR$lossOrder),  # 16
			as.integer(0),			  # 17
			PACKAGE="ppMeasures")
		tR$distance <- cOut[[13]]
		tR$xtoy <- cOut[[10]]+1
		tR$ytox <- cOut[[11]]+1
		tR$algorithm <- 'IMA'
	}
	return(tR)
}


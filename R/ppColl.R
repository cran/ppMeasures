ppColl <-
function(pts, key, wts=NULL, wtsKey=NULL, nMissing=0){
	tR <- list()
	class(tR) <- 'ppColl'
	tR$points <- NA
	tR$pattID <- key
	tR$key    <- NA
	tR$keyMax <- NA
	tR$wts    <- NA
	tR$dim    <- NA
	tR$maxObs <- NA
	tR$obs    <- NA
	
	if(is.data.frame(pts)){
		pts <- as.matrix(pts)
	}
	if(is.vector(pts)){
		tR$dim    <- 1
		tR$points <- matrix(pts, ncol=1)
	} else {
		tR$points <- pts
		tR$dim    <- dim(tR$points)[2]
	}
	
	if(!is.null(wts[1])){
		if(is.null(wtsKey[1])){
			stop('Must provide "wtsKey" if providing "wts".\n')
		}
		if(length(wts) != length(wtsKey)){
			stop('Arguments "wts" and "wtsKey" must have same length.\n')
		}
		if(!all(key %in% wtsKey)){
			stop('Weights are missing for some patterns.\n')
		}
		tR$key <- match(tR$pattID, wtsKey)
		tR$keyMax <- length(wtsKey)
		tR$wts <- wts
	} else {
		keyWts <- unique(tR$pattID)
		tR$key  <- match(tR$pattID, keyWts)
		tR$keyMax <- max(tR$key)+nMissing
		tR$wts <- rep(1, tR$keyMax)
	}
	
	tR$maxObs <- 0
	for(i in 1:tR$keyMax){
		tR$obs[i] <- sum(tR$key == i)
		if(tR$obs[i] > tR$maxObs){
			tR$maxObs <- tR$obs[i]
		}
	}
	
	if(length(tR$key) == dim(tR$points)[1]){
		return(tR)
	} else {
		stop('The number of points and keys must match.\n')
	}
}


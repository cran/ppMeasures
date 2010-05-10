points.ppColl <-
function(x, pch='default', cex='propToWts', col='default', Dims=1:2, ...){
	if(x$dim > 2){
		hold     <- x
		if(length(Dims) > 1){
			hold$points <- x$points[,Dims[1:2]]
			hold$dim    <- 2
		} else {
			hold$points <- matrix(x$points[,Dims[1]], ncol=1)
			hold$dim    <- 1
		}
		x <- hold
	}
	if(pch[1] == 'default'){
		if(x$dim == 1){
			pch <- 20
		} else {
			pch <- x$key
			pch <- ((pch-1) %% 25)+1
		}
	} else if(length(pch) == x$keyMax){
		pch <- pch[x$key]
	} else {
		pch <- rep(pch, length(x$key))
	}
	if(col[1] == 'default'){
		if(x$dim == 1){
			col <- '#00000088'
		} else {
			col <- x$key
		}
	} else if(length(col) == x$keyMax){
		col <- col[x$key]
	} else {
		col <- rep(col, length(x$key))
	}
	if(cex[1] == 'propToWts'){
		cex <- sqrt(x$wts[x$key]/min(x$wts[x$key])*0.7)
	} else if(length(cex) < x$keyMax){
		cex <- rep(cex, length(x$key))
	} else if(length(cex) == x$keyMax){
		cex <- cex[x$key]
	}
	if(x$dim == 1 || length(Dims) == 1){
		if(x$dim == 1){
			Dims <- 1
		}
		points(x$points[,Dims], x$key, pch=pch, cex=cex, col=col, ...)
	} else {
		points(x$points, pch=pch, cex=cex, col=col, ...)
	}
	nLeft <- x$keyMax - length(unique(x$key))
	if(nLeft > 0){
		cat('Number of patterns with no points:',
			nLeft, '\n')
	}
}


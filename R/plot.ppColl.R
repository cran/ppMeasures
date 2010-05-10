plot.ppColl <-
function(x, pch='default', cex='propToWts', col='default', xlab=NULL, ylab=NULL, Dims=1:2, ylim=NULL, addLines=FALSE, ...){
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
	if(x$dim == 1){
		if(is.null(xlab)[1]){
			xlab <- paste('dimension', Dims[1])
		}
		if(is.null(ylab)[1]){
			ylab <- 'pattern'
		}
		if(is.null(ylim[1])){
			R <- range(x$key)
			ylim <- R + c(-1,1)*diff(R)/10
		}
		plot(x$points[,1], x$key, pch=pch, cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
		if(addLines[1]){
			abline(h=0:x$keyMax, col='#00000033')
			abline(h=seq(0,x$keyMax,5), col='#00000066')
		}
	} else {
		if(is.null(xlab)[1]){
			xlab <- paste('dim', Dims[1])
		}
		if(is.null(ylab)[1]){
			ylab <- paste('dim', Dims[2])
		}
		plot(x$points, pch=pch, cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
	}
	nLeft <- x$keyMax - length(unique(x$key))
	if(nLeft > 0){
		cat('Number of patterns with no points:',
			nLeft, '\n')
	}
}

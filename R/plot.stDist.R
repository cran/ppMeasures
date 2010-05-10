plot.stDist <-
function(x, Dims=1:2, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, pch=c(1,20), col=c('#00000088', '#FF000088'), cex=c(1,1.5), yPos=0:1, legend=FALSE, legPos='topright', ...){
	if(x$dim > 2){
		hold     <- x
		if(length(Dims) > 1){
			hold$x   <- x$x[,Dims[1:2]]
			hold$y   <- x$y[,Dims[1:2]]
			hold$dim <- 2
		} else {
			hold$x   <- matrix(x$x[,Dims[1]], ncol=1)
			hold$y   <- matrix(x$y[,Dims[1]], ncol=1)
			hold$dim <- 1
		}
		x <- hold
	}
	if(length(pch) == 1){
		pch <- rep(pch, 2)
	}
	if(length(col) == 1){
		col <- rep(col, 2)
	}
	if(length(cex) == 1){
		cex <- rep(cex, 2)
	}
	if(x$algorithm == 'MSU'){
		stop('Cannot plot an stDist object when using the MSU algorithm.\n')
	}
	if(x$dim == 2){
		if(is.null(xlim)[1]){
			xlim <- range(c(x$x[,1], x$y[,1]))
		}
		if(is.null(ylim)[1]){
			ylim <- range(c(x$x[,2], x$y[,2]))
		}
		if(is.null(xlab)){
			xlab <- paste('dim', Dims[1])
		}
		if(is.null(ylab)){
			ylab <- paste('dim', Dims[2])
		}
		plot(x$x[,1:2], col=col[1], pch=pch[1], cex=cex[1],
			xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
		points(x$y[,1:2], col=col[2], pch=pch[2], cex=cex[2])
		for(i in 1:length(x$xtoy)){
			if(x$xtoy[i] > 0){
				lines(c(x$x[i,1], x$y[x$xtoy[i],1]),						c(x$x[i,2], x$y[x$xtoy[i],2]))
			}
		}
	} else {
		if(is.null(xlab)){
			xlab <- paste('dimension', Dims[1])
		}
		if(is.null(ylab)){
			ylab <- 'pattern'
		}
		if(is.null(xlim)[1]){
			xlim <- range(c(x$x[,1], x$y[,1]))
		}
		yx <- rep(yPos[1], dim(x$x)[1])
		yy <- rep(yPos[2], dim(x$y)[1])
		if(is.null(ylim)[1]){
			ylim <- sort(yPos)+abs(diff(yPos))*0.2*c(-1,1)
		}
		plot(x$x[,1], yx, col=col[1], pch=pch[1],
			cex=cex[1], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
		points(x$y[,1], yy, col=col[2], pch=pch[2],
			cex=cex[2])
		for(i in 1:length(x$xtoy)){
			if(x$xtoy[i] > 0){
				lines(c(x$x[i], x$y[x$xtoy[i]]), yPos)
			}
		}
	}
	if(legend[1] != FALSE){
		legend(legPos, pch=pch, col=col, pt.cex=cex, legend=legend)
	}
}


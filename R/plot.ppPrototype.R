plot.ppPrototype <-
function(x, xlim=NULL, ylim=NULL, pch=1, col=1, cex=1, yPos=NULL, Dims=1:2, xlab=NULL, ylab=NULL, ...){
	if(x$dim >= 2){
		hold           <- x
		hold$prototype <- x$prototype[,Dims[1:2]]
		hold$dim       <- 2
		x              <- hold
		if(is.null(xlab)[1]){
			xlab <- paste('dim', Dims[1])
		}
		if(is.null(ylab)[1]){
			ylab <- paste('dim', Dims[2])
		}
		if(is.null(xlim)[1]){
			xlim <- range(x$prototype[,1])
		}
		if(is.null(ylim)[1]){
			ylim <- range(x$prototype[,2])
		}
		plot(x$prototype[,1], x$prototype[,2], col=col, pch=pch, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	} else {
		if(is.null(xlab)[1]){
			xlab <- paste('dimension', Dims[1])
		}
		if(is.null(ylab)[1]){
			ylab <- ''
		}
		if(is.null(xlim)[1]){
			xlim <- range(x$prototype[,1])
		}
		if(is.null(yPos[1])){
			yPos <- rep(0, dim(x$prototype)[1])
		} else if(length(yPos) == 1){
			yPos <- rep(yPos, dim(x$prototype)[1])
		}
		if(is.null(ylim)[1]){
			ylim <- yPos[1] + 0.5*c(-1,1)
		}
		plot(x$prototype[,Dims[1]], yPos, col=col, pch=pch, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	}
}


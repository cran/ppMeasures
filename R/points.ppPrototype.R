points.ppPrototype <-
function(x, Dims=1:2, at=0, ...){
	if(x$dim == 1){
		points(x$prototype, rep(at, length(x$prototype)), ...)
	} else {
		hold           <- x
		hold$prototype <- x$prototype[,Dims[1:2]]
		hold$dim       <- 2
		x              <- hold
		points(x$prototype[,1], x$prototype[,2], ...)
	}
}


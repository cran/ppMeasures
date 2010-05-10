print.ppColl <-
function(x, ...){
	pts <- x$points
	pts <- cbind(x$key, pts)
	colnames(pts) <- c('key', paste('dim', 1:x$dim))
	rownames(pts) <- 1:dim(pts)[1]
	print(pts)
}


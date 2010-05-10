summary.stDist <-
function(object, ...){
	nx <- length(object$x)
	ny <- length(object$y)
	n  <- sum(object$xtoy > 0)
	cat('Algorithm: ', object$algorithm, '\n', sep='')
	if(object$algorithm == 'IMA'){
		cat('Max branch:', object$maxBranch, '\n')
	}
	if(object$algorithm != 'MSU'){
		if(n != 1){
			cat(n, 'points were matched\n')
		} else {
			cat(n, 'point was matched\n')
		}
	}
	cat('Distance:', object$distance, '\n')
}


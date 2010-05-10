summary.ppColl <-
function(object, ...){
	n <- min(c(length(object$points), dim(object$points)[1]))
	if(n != 1){
		cat('There are', n, 'points over ')
	} else {
		cat('There is', n, 'point over ')
	}
	if(object$dim > 1){
		cat(object$dim, 'dimensions and ')
	} else {
		cat(object$dim, 'dimension and ')
	}
	if(object$keyMax > 1){
		cat(object$keyMax, 'patterns.\n')
	} else {
		cat(object$keyMax, 'pattern.\n')
	}
	nMissing <- object$keyMax-length(unique(object$key))
	if(nMissing > 0){
		cat('Number of patterns with no points: ')
		cat(nMissing, '.\n', sep='')
	} else {
		cat('All patterns contain at least one point.\n')
	}
	if(all(object$wts == object$wts[1])){
		cat('All patterns are weighted equally.\n\n')
	} else {
		cat('Patterns have different weights.\n\n')
	}
	
	R <- apply(object$points, 2, fivenum)
	colnames(R) <- paste('dim', 1:object$dim)
	rownames(R) <- c('min', 'Q1', 'median', 'Q3', 'max')
	print(signif(R, 4))
}


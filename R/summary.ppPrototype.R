summary.ppPrototype <-
function(object, ...){
	cat('A prototype was computed based on', object$ppColl$keyMax, 'patterns.\n')
	cat('pm: ')
	cat(object$pm, sep=', ')
	cat('\n')
	cat('pa:', object$pa, '\n')
	cat('pd:', object$pd, '\n')
	cat('Total cost:', object$cost, '\n\n')
	
	toPrint <- object$prototype
	colnames(toPrint) <- paste('dim', 1:dim(toPrint)[2])
	rownames(toPrint) <- 1:dim(toPrint)[1]
	print(signif(toPrint, 4))
}


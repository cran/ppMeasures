print.ppPrototype <-
function(x, ...){
	if(!is.null(x$prototype)[1]){
		toPrint <- matrix(x$prototype,
			length(x$prototype)/x$dim, x$dim)
		colnames(toPrint) <- paste('dim', 1:dim(toPrint)[2])
		rownames(toPrint) <- 1:dim(toPrint)[1]
		print(signif(toPrint, 4))
	} else {
		cat("There are no prototype points.")
	}
}



#################################################################
# transforming data in BIFIE.data object
BIFIE.data.transform <- function( bifieobj , transform.formula ){
	datalistM <- bifieobj$datalistM
	varnames <- bifieobj$varnames
	dfr <- as.data.frame(datalistM)
	colnames(dfr) <- varnames
	M1 <- model.matrix( transform.formula , dfr )
	varnames.added <- colnames(M1)
	N1 <- nrow(datalistM)
	N2 <- ncol(datalistM)
	varsindex.added <- seq( N2 + 1 , N2 + ncol(M1) )
	M1.new <- matrix( NA , nrow=N1 , ncol=ncol(M1) )
	colnames(M1.new) <- varnames.added
	VV <- ncol(M1)
	M1.new[ match( rownames(M1),rownames(dfr) ) , ] <- M1
	datalistM <- as.matrix( cbind( datalistM , M1.new ) )
	colnames(datalistM) <- NULL
	varnames1 <- c( varnames , varnames.added )
	#****
	# dat1
	N <- bifieobj$N
	Nimp <- bifieobj$Nimp
	dat1 <- datalistM[ seq( N*(Nimp-1) + 1 , Nimp*N ) , ]
	colnames(dat1) <- varnames1
	bifieobj$dat1 <- as.matrix( dat1)
	bifieobj$datalistM <- datalistM
	bifieobj$varnames <- varnames1
	bifieobj$varnames.added <- varnames.added
	bifieobj$varsindex.added <- varsindex.added
	cat( paste0( "Included " , VV , " variables: " , paste0( varnames.added , collapse=" " ) , "\n") )
	return( bifieobj )
		}
#################################################################		
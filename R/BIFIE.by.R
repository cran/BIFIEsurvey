

#######################################################################
# BIFIE.by function
BIFIE.by <- function( BIFIEobj , vars , userfct , userparnames=NULL ,
		group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	bifieobj <- BIFIEobj	
	if (bifieobj$cdata){
		varnames <- unique( c( vars , group , "one") )
		bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames=varnames )	
						}				
	FF <- Nimp <- bifieobj$Nimp
	N <- bifieobj$N
	dat1 <- bifieobj$dat1
	wgt <- bifieobj$wgt
	wgtrep <- bifieobj$wgtrep
	varnames <- bifieobj$varnames
	RR <- bifieobj$RR
	datalistM <- bifieobj$datalistM
    fayfac <- bifieobj$fayfac	
	
	if (RR==1){ RR <- 0 }
	if ( ! se ){ 
		wgtrep <- matrix( wgt , ncol=1 )
		RR <- 0
				}	
	
	vars_index <- unlist( sapply( vars , FUN = function(vv){ 
						which( varnames == vv ) } ) )
    # vars values
	VV <- length(vars)
					
	wgt_ <- matrix( wgt , ncol=1 )
	if ( is.null( group) ){ nogroup <- TRUE } else { nogroup <- FALSE }
	cat(paste0( "|" , paste0( rep("*" , FF) , collapse="") , "|\n" ))
	if (nogroup){
	    group <- "one"
	    group_values <- c(1)
			}
    group_index <- which( varnames %in% group )
    if ( is.null(group_values ) ){ 
		t1 <- table( dat1[ , group_index ] )				  
	    group_values <- sort( as.numeric( paste( names(t1) ) ))
				}
				
	#**************************************************************************#
	# Rcpp call

	res <- .Call("bifie_by" , datalistM , wgt_ , wgtrep ,	vars_index - 1,    fayfac ,
				Nimp , group_index - 1 , group_values , userfct , PACKAGE="BIFIEsurvey")
	NP <- res$NP
	GG <- length(group_values)
	ZZ <- NP
	if (is.null( userparnames ) ){
		userparnames <- paste0("parm",1:NP) 
				}
	
	dfr <- data.frame( "parm" = rep( userparnames , GG )
						)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , each=ZZ )
	             }				 				 	



	dfr$Ncases <- rep( rowMeans( res$ncasesM ) , each=ZZ )
	dfr$Nweight <- rep( rowMeans( res$sumwgtM ) , each=ZZ )	
	dfr$est <- res$parsL$pars
	dfr$SE <- res$parsL$pars_se
	dfr$fmi <- res$parsL$pars_fmi
	dfr$VarMI <- res$parsL$pars_varBetween
	dfr$VarRep <- res$parsL$pars_varWithin
	if ( ( ! se ) &  ( RR==0 ) ){				
		dfr$SE <- dfr$fmi <- dfr$VarMI <- dfr$VarRep <- NULL
				}				


	# create vector of parameter names
	parnames <- paste0( dfr$parm   , "_" , dfr$group , dfr$groupval )

	
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr , 
			"output" = res , 	"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac , "GG"=GG ,
			"parnames" = parnames)
	class(res1) <- "BIFIE.by"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.by function

summary.BIFIE.by <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Statistical Inference for User Definition Function \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}
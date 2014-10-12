

#######################################################################
# Correlations and covariances
BIFIE.correl <- function( BIFIEobj , vars , group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	cl19 <- match.call()
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
		t1 <- fasttable( datalistM[ , group_index ] )				  
	    group_values <- sort( as.numeric( paste( names(t1) ) ))
				}
				
	#**************************************************************************#
	# Rcpp call

	res <- .Call( "bifie_correl" , datalistM , wgt_ , as.matrix(wgtrep) , vars_index -1 , fayfac , 
                Nimp ,  group_index -  1, group_values , PACKAGE="BIFIEsurvey" )
	GG <- length(group_values)
    itempair_index <- res$itempair_index + 1
	ZZ <- nrow(itempair_index )	

	dfr <- data.frame( "var1" = rep( vars[ itempair_index[,1] ] , each=GG ) ,
					"var2" = rep( vars[ itempair_index[,2] ] , each=GG )
						)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , ZZ )
	             }				 				 				 	 
	dfr$Ncases <- rep( rowMeans( res$ncases1M ) , ZZ )
	dfr$Nweight <- rep( rowMeans( res$sumwgt1M ) , ZZ )	
	dfr$cor <- res$cor1$pars
	dfr$cor_SE <- res$cor1$pars_se
	dfr$cor_fmi <- res$cor1$pars_fmi
	dfr$cor_VarMI <- res$cor1$pars_varBetween
	dfr$cor_VarRep <- res$cor1$pars_varWithin
	if (RR==0){				
		dfr$cor_SE <- dfr$cor_fmi <- dfr$cor_VarMI <- dfr$cor_VarRep <- NULL
				}
	#	i1 <- match( dfr$var1 , vars )
	#	i2 <- match( dfr$var2 , vars )	
	dfr <- dfr[ dfr$var1 != dfr$var2 , ]	
	dfr.cor <- dfr			
				

	dfr <- data.frame( "var1" = rep( vars[ itempair_index[,1] ] , each=GG ) ,
					"var2" = rep( vars[ itempair_index[,2] ] , each=GG )
						)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , ZZ )
	             }				 				 				 	 
	dfr$Ncases <- rep( rowMeans( res$ncases1M ) , ZZ )
	dfr$Nweight <- rep( rowMeans( res$sumwgt1M ) , ZZ )	
	dfr$cov <- res$cov1$pars
	dfr$cov_SE <- res$cov1$pars_se
	dfr$cov_fmi <- res$cov1$pars_fmi
	dfr$cov_VarMI <- res$cov1$pars_varBetween
	dfr$cov_VarRep <- res$cov1$pars_varWithin
	if ( ( ! se ) &  ( RR==0 ) ){				
		dfr$cov_SE <- dfr$cov_fmi <- dfr$cov_VarMI <- dfr$cov_VarRep <- NULL
				}				
	if ( Nimp == 1){				
		dfr$cov_fmi <- dfr$cov_VarMI <- NULL
				}					
				
	dfr.cov <- dfr							

	#*****
	# construct list of correlation matrices
	ml <- as.list(1:GG)
	names(ml) <- paste0(group,group_values)
	#*** correlation matrix
	ml0 <- ml
	cl <- res$cor1_matrix
	for (gg in 1:GG){
		ml0[[gg]] <- cl[ , 1:VV + (gg-1 )*VV ]
		colnames(ml0[[gg]]) <- rownames(ml0[[gg]]) <- vars 
					}		
	cor_matrix <- ml0
	#*** covariance matrix
	ml0 <- ml
	cl <- res$cov1_matrix
	for (gg in 1:GG){
		ml0[[gg]] <- cl[ , 1:VV + (gg-1 )*VV ]
		colnames(ml0[[gg]]) <- rownames(ml0[[gg]]) <- vars 
					}		
	cov_matrix <- ml0	

	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )		
	parnames <- paste0( dfr$var1   , "_" , dfr$var2 , 			
			ifelse( ! nogroupL , paste0("_" , dfr$groupvar ) , "" ) ,
			ifelse( ! nogroupL , paste0( "_" , dfr$groupval  ) , "" ) 
			)
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.cor" = dfr.cor , "stat.cov" = dfr.cov , 
			"output" = res , "cor_matrix" =  cor_matrix , 
			"cov_matrix" = cov_matrix , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"itempair_index" = itempair_index , "GG"=GG ,
			"parnames" = parnames , "CALL"= cl19)
	class(res1) <- "BIFIE.correl"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.correl function
summary.BIFIE.correl <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Statistical Inference for Correlations \n")	
	obji <- object$stat.cor
	print.object.summary( obji , digits=digits )
	cat("\nCorrelation Matrices \n\n")	
	obji <- object$cor_matrix
	GG <- object$GG
	for (gg in 1:GG){
	    obji[[gg]] <- round( obji[[gg]] , digits=digits)
				}
	print(obji)
#	print.object.summary( obji , digits=digits )			
			}


#######################################################################
# frequency tables
BIFIE.freq <- function( BIFIEobj , vars , group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	cl <- match.call()		
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
	vars_info <- list(1:VV)
	for (vv in 1:VV){
	   # t1 <- table( dat1[,vars_index[vv] ] )	   
    	t1 <- fasttable( datalistM[ , vars_index[vv] ] )
	   vars_info[[vv]] <- sort( as.numeric( paste0(names(t1) )))	   
		    }
	vars_values_numb <- unlist( lapply( vars_info , FUN = function(uu){ length(uu) } )	) 
	vars_values <- matrix(NA, nrow=max(vars_values_numb) , ncol=VV)
	for (vv in 1:VV){
	   vars_values[ seq(1,vars_values_numb[vv] ) , vv ] <- vars_info[[vv]]
			}
					
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

	res <-  .Call("bifie_freq"  ,datalistM , wgt_ , as.matrix(wgtrep) , vars_index -1 , fayfac , 
				Nimp ,  group_index -  1, group_values , as.matrix(vars_values) ,
				vars_values_numb , PACKAGE="BIFIEsurvey" )

	GG <- res$outlist$GG

	dfr <- data.frame( "var" = rep( rep( vars , vars_values_numb ) , each=GG ) )
	VV <- length(vars)
	varval <- unlist( sapply( 1:VV , FUN = function(vv){
		# vv <- 1
		rep( vars_values[ 1:vars_values_numb[vv] , vv ] , GG )
				} ) )
	dfr$varval <- varval
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( rep( group_values , VV) , rep(vars_values_numb,each=GG) )
	             }
	dfr$Ncases <- rowMeans( res$ncases1M )
	dfr$Nweight <- res$perc1$pars
	# percentage
	dfr$perc <- res$perc2$pars
	dfr$perc_SE <- res$perc2$pars_se
	dfr$perc_fmi <- res$perc2$pars_fmi
	dfr$perc_VarMI <- res$perc2$pars_varBetween
	dfr$perc_VarRep <- res$perc2$pars_varWithin
			

	if ( ( ! se ) &  ( RR==0 ) ){				
		dfr$perc_SE <- dfr$perc_fmi <- dfr$perc_VarMI <- dfr$perc_VarRep <- NULL
				}				
	if ( Nimp==1 ){				
		dfr$perc_fmi <- dfr$perc_VarMI <- NULL
				}		
	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )
	parnames <- paste0( dfr$var   , "_" , dfr$varval , 
			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
			ifelse( ! nogroupL , dfr$groupval , "" ) )	
	
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr , "output" = res , "timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"parnames" = parnames , "CALL"=cl )
	class(res1) <- "BIFIE.freq"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.freq function

summary.BIFIE.freq <- function( object , digits=3 , ... ){
    BIFIE.summary(object)
	cat("Relative Frequencies \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}
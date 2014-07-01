
#######################################################################
# univariate statistics
BIFIE.univar <- function( BIFIEobj , vars , group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	bifieobj <- BIFIEobj
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
	#****************** no grouping variable **********************************#
	if ( nogroup ){
#		res <- .Call( "univar_1group" ,  datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac ,  Nimp ,
#				PACKAGE="BIFIEsurvey" ) 
#		GG <- 1
#		VV <- length(vars)
		res <- .Call( "univar_multiple_V2group" , datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac , Nimp ,
				group_index - 1, group_values , PACKAGE="BIFIEsurvey" )								
		GG <- length(group_values)
		VV <- length(vars)	
		dfr <- data.frame( "var" = rep(vars,each=GG) , 
				"Nweight" = rowMeans(res$sumweightM) ,
				"Ncases" = rowMeans( res$ncasesM) , 
				"M" = res$mean1 , "M_SE" = res$mean1_se ,
				"M_fmi" = res$mean1_fmi   , "M_VarMI"= res$mean1_varBetween , "M_VarRep"= res$mean1_varWithin   ,
				"SD" = res$sd1 , "SD_SE" = res$sd1_se ,
				"SD_fmi" = res$sd1_fmi   , "SD_VarMI"= res$sd1_varBetween , "SD_VarRep"= res$sd1_varWithin            
					)
				}
	#**************************************************************************#
	#****************** with grouping variable ********************************#		
	if ( ! nogroup ){
		res <- .Call( "univar_multiple_V2group" , datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac , Nimp ,
				group_index - 1, group_values , PACKAGE="BIFIEsurvey" )								
		GG <- length(group_values)
		VV <- length(vars)
		dfr <- data.frame( "var" = rep(vars,each=GG) , 
				"groupvar" = group , 
				"groupval" = rep(group_values  , VV ) , 
				"Nweight" = rep( rowMeans(res$sumweightM) , VV ) ,
				"Ncases" = res$ncases ,
				"M" = res$mean1 , "M_SE" = res$mean1_se ,
				"M_fmi" = res$mean1_fmi   , "M_VarMI"= res$mean1_varBetween , "M_VarRep"= res$mean1_varWithin   ,
				"SD" = res$sd1 , "SD_SE" = res$sd1_se ,
				"SD_fmi" = res$sd1_fmi   , "SD_VarMI"= res$sd1_varBetween , "SD_VarRep"= res$sd1_varWithin            
					 )
	         	}

	if (RR==0){				
		dfr$M_SE <- dfr$M_fmi <- dfr$M_VarMI <- dfr$M_VarRep <- NULL
		dfr$SD_SE <- dfr$SD_fmi <- dfr$SD_VarMI <- dfr$SD_VarRep <- NULL		
				}				
	
	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )
	parnames <- paste0( dfr$var   , 
			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
			ifelse( ! nogroupL , dfr$groupval , "" ) )		
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr , "output" = res , "timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac , "parnames" = parnames ,
			"GG" = GG , "VV"=VV , "vars" = vars , "group" = group )
	class(res1) <- "BIFIE.univar"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.univar function

summary.BIFIE.univar <- function( object , digits=3 , ... ){
    BIFIE.summary(object)     
	cat("Univariate Statistics\n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}
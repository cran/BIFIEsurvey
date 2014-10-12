

#######################################################################
# BIFIE Wald test
BIFIE.waldtest <- function( BIFIE.method , Cdes , rdes , type=NULL ){
	#****
	s1 <- Sys.time()
	cl <- match.call()
	res1 <- BIFIE.method

	# extract replicated parameters
	parsres <- extract.replicated.pars( BIFIE.method = res1 , type=type)	
	parsM <- parsres$parsM	
	parsrepM <- parsres$parsrepM	
	
	fayfac <- res1$fayfac
	
	#****************************************			
	# which columns in C do have non-zero entries
	Ccols <- which( colSums( abs( Cdes) ) > 0 )

	# apply Rcpp Wald test function
	if (TRUE){
	   res <- .Call("bifie_waldtest" ,  parsM , parsrepM , Cdes , rdes , Ccols - 1 , fayfac ,
			           PACKAGE="BIFIEsurvey")					   
				}
#    if (FALSE){
#	  res <- bifie_waldtest(  parsM , parsrepM , Cdes , rdes , Ccols - 1 , fayfac )				
#	          }
	RR <- res$RR
    Nimp <- res$Nimp
	fayfac <- res$fayfac
	N <- BIFIE.method$N
	# data frame with results
	dfr <- data.frame( "D1"= res$D1 , "D2" = res$D2 , "df1"= res$df , 
			"D1_df2"= round(res$nu2,1) , "D2_df2" = round(res$nu3,1) ,
			"D1_p" = res$p_D1 , "D2_p" = res$p_D2 )
	# rownames(dfr) <- NULL			
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.D" = dfr ,
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"class.BIFIE.method" = class(BIFIE.method) , "CALL"= cl )
	class(res1) <- "BIFIE.waldtest"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.correl function
summary.BIFIE.waldtest <- function( object , digits=4 , ... ){
    BIFIE.summary(object , FALSE)
	cat("D1 and D2 Statistic for Wald Test \n\n")	
	obji <- object$stat.D
	print.object.summary( obji , digits=digits )			
			}
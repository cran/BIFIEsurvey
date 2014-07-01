

#######################################################################
# BIFIE Wald test
BIFIE.waldtest <- function( BIFIE.method , Cdes , rdes ){
	#****
	s1 <- Sys.time()
	res1 <- BIFIE.method
	#************************************
	#**** linear regression
	if ( class( BIFIE.method) == "BIFIE.linreg"){ 	
		# parameters in every imputed dataset
		parsM <- res1$output$regrcoefM
		# replicated parameters
		parsrepM <- res1$output$regrcoefrepM
				}
	#************************************
	#**** correlation
	if ( class( BIFIE.method) == "BIFIE.correl"){ 	
		parsM <- res1$output$cor1M
		parsrepM <- res1$output$cor1repM
				}				
	#************************************
	#**** frequencies
	if ( class( BIFIE.method) == "BIFIE.freq"){ 	
		parsM <- res1$output$perc2M
		parsrepM <- res1$output$perc2repM
				}		
	#************************************
	#**** univar
	if ( class( BIFIE.method) == "BIFIE.univar"){ 	
		parsM <- res1$output$mean1M
		parsrepM <- res1$output$mean1repM
				}		
	#************************************
	#**** crosstab
	if ( class( BIFIE.method) == "BIFIE.crosstab"){ 	
		parsM <- res1$output$ctparsM
		parsrepM <- res1$output$ctparsrepM
				}

				
	fayfac <- res1$fayfac
	
	#****************************************			
	# which columns in C do have non-zero entries
	Ccols <- which( colSums( abs( Cdes) ) > 0 )

	# apply Rcpp Wald test function
	if (TRUE){
	   res <- .Call("bifie_waldtest" ,  parsM , parsrepM , Cdes , rdes , Ccols - 1 , fayfac ,
			           PACKAGE="BIFIEsurvey")
				}
    if (FALSE){
	  res <- bifie_waldtest(  parsM , parsrepM , Cdes , rdes , Ccols - 1 , fayfac )				
	          }
	# ariv <- sum( diag( res$var_b %*% solve( res$var_w ) ) )
	# ariv <- ariv * ( 1 + 1/res$Nimp ) / res$df
	# var_t <- ( 1 + ariv) * res$var_w
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
			"class.BIFIE.method" = class(BIFIE.method) )
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
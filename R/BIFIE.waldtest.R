

#######################################################################
# BIFIE Wald test
BIFIE.waldtest <- function( BIFIE.method , Cdes , rdes , type=NULL ){
	#****
	s1 <- base::Sys.time()
	cl <- base::match.call()
	res1 <- BIFIE.method

	# extract replicated parameters
	parsres <- extract.replicated.pars( BIFIE.method = res1 , type=type)	
	parsM <- parsres$parsM	
	parsrepM <- parsres$parsrepM
    parnames <- parsres$parnames	
	fayfac <- res1$fayfac
	N <- BIFIE.method$N	
	Nimp <- BIFIE.method$Nimp
	RR <- BIFIE.method$RR
	
	#****************************************			
	# which columns in C do have non-zero entries
	Ccols <- base::which( base::colSums( base::abs( Cdes) ) > 0 )

	if ( ! BIFIE.method$NMI ){
		
		# apply Rcpp Wald test function
		if (TRUE){
		   res <- base::.Call("bifie_waldtest" ,  parsM   , parsrepM , Cdes , rdes , Ccols - 1 , fayfac ,
						   PACKAGE="BIFIEsurvey")					   
					}
					
	#    if (FALSE){
	#	  res <- bifie_waldtest(  parsM , parsrepM , Cdes , rdes , Ccols - 1 , fayfac )				
	#	          }
		RR <- res$RR
		Nimp <- res$Nimp
		fayfac <- res$fayfac
		
		# data frame with results
		dfr <- base::data.frame( "D1"= res$D1 , "D2" = res$D2 , "df1"= res$df , 
				"D1_df2"= base::round(res$nu2,1) , "D2_df2" = base::round(res$nu3,1) ,
				"D1_p" = res$p_D1 , "D2_p" = res$p_D2 )
			}
			
	if ( BIFIE.method$NMI ){
		Cdes_cols <- Cdes[ , Ccols , drop=FALSE]
		df1 <- nrow(Cdes_cols)
		parsM2 <- Cdes_cols %*% parsM[ Ccols , ] 
		parsrepM2 <- Cdes_cols %*% parsrepM[ Ccols , ]
		# within covariance matrices
		res0 <- base::.Call( "bifie_comp_vcov_within" , parsM2 , parsrepM2 , fayfac , 
					BIFIE.method$RR , Nimp , package="BIFIEsurvey" )
		u <- res0$u
		Nimp_NMI <- BIFIE.method$Nimp_NMI
		qhat <- base::array( parsM2 , dim= c( df1 , Nimp_NMI[2] , Nimp_NMI[1] ) )
		qhat <- base::aperm( qhat , c(3,2,1) )
		v1 <- base::paste0("parm",1:df1)
		dimnames(qhat) <- base::list( 
			paste0("imp_nmi_dim1_" , seq(1,dim(qhat)[[1]] ) ) ,
			paste0("imp_nmi_dim2_" , seq(1,dim(qhat)[[2]] ) ) ,
			v1 )
	
		if ( ! is.null( dimnames(qhat) ) ){ 
			dimnames(qhat)[[3]] <- v1
											}
		u <- base::array( u , dim = c( df1 , df1 , Nimp_NMI[2] , Nimp_NMI[1] ) )
		u <- base::aperm( u , c(4,3,1,2) )	
		res <- miceadds::NMIwaldtest( qhat= qhat , u = u , testnull=v1)
		# res <- NMIwaldtest( qhat= qhat , u = u , testnull=v1)
		dfr <- base::data.frame( "D1"= res$stat$F ,  "df1"= res$stat$df1 , 
				"D1_df2"= base::round(res$stat$df2,1) , 
				"D1_p" = res$stat$pval  )
				}
		
	# rownames(dfr) <- NULL			
	
	#*************************** OUTPUT ***************************************
	s2 <- base::Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- base::list( "stat.D" = dfr ,
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIE.method$NMI , "Nimp_NMI" = BIFIE.method$Nimp_NMI , 
			"class.BIFIE.method" = class(BIFIE.method) , "CALL"= cl )
	base::class(res1) <- "BIFIE.waldtest"
	base::return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.correl function
summary.BIFIE.waldtest <- function( object , digits=4 , ... ){
    BIFIE.summary(object , FALSE)
	if ( ! object$NMI ){ base::cat("D1 and D2 Statistic for Wald Test \n\n")	 }
	if (  object$NMI ){ base::cat("D1 Statistic for Wald Test \n\n")	 }
	obji <- object$stat.D
	print.object.summary( obji , digits=digits )			
			}
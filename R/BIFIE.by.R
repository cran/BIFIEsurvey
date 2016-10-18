

#######################################################################
# BIFIE.by function
BIFIE.by <- function( BIFIEobj , vars , userfct , userparnames=NULL ,
		group=NULL , group_values=NULL , se=TRUE , use_Rcpp = TRUE ){
	#****
	s1 <- base::Sys.time()
	cl <- base::match.call()
	bifieobj <- BIFIEobj	
	if (bifieobj$cdata){
		varnames <- base::unique( c( vars , group , "one") )
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
		wgtrep <- base::matrix( wgt , ncol=1 )
		RR <- 0
				}	
	
	vars_index <- base::unlist( base::sapply( vars , FUN = function(vv){ 
						base::which( varnames == vv ) } , simplify=TRUE ) )
    # vars values
	VV <- base::length(vars)
					
	wgt_ <- base::matrix( wgt , ncol=1 )
	if ( is.null( group) ){ nogroup <- TRUE } else { nogroup <- FALSE }
	base::cat( base::paste0( "|" , base::paste0( base::rep("*" , FF) , collapse="") , "|\n" ))
	if (nogroup){
	    group <- "one"
	    group_values <- c(1)
			}
			

	#@@@@***
    group_index <- base::match( group , varnames )
	#@@@@***

    if ( is.null(group_values ) ){ 
		t1 <- fasttable( datalistM[ , group_index ] )				  
	    group_values <- base::sort( base::as.numeric( base::paste( base::names(t1) ) ))
				}
	
	#@@@@***
	res00 <- BIFIE_create_pseudogroup( datalistM , group , group_index , group_values )				
	res00$datalistM -> datalistM 
	res00$group_index -> group_index
	res00$GR -> GR 
	res00$group_values -> group_values
	res00$group -> group
	#@@@@***			

		
	#****
	# pure R implementation
	if ( ! use_Rcpp ){	
		res <- BIFIE_by_helper_pureR(
			group_values , userfct , datalistM ,
			N , vars_index , wgt_ , wgtrep , Nimp , RR , fayfac ,
			group_index , userparnames
				)
					}
	
	#****
	# Rcpp implementation
	if ( use_Rcpp ){
		res <- base::.Call("bifie_by" , datalistM , wgt_ , wgtrep ,	vars_index - 1,    fayfac ,
				Nimp , group_index - 1 , group_values , userfct , PACKAGE="BIFIEsurvey")
					}
	
	NP <- res$NP
	GG <- base::length(group_values)
	ZZ <- NP
	if (is.null( userparnames ) ){
		userparnames <- base::paste0("parm",1:NP) 
				}
	
	dfr <- base::data.frame( "parm" = base::rep( userparnames , GG )
							)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- base::rep( group_values , each=ZZ )
	             }				 				 	


	dfr$Ncases <- base::rep( rowMeans( res$ncasesM ) , each=ZZ )
	dfr$Nweight <- base::rep( rowMeans( res$sumwgtM ) , each=ZZ )

	dfr <- create_summary_table( res_pars=res$parsL , 
				     parsM=res$parsM   , parsrepM=res$parsrepM , 
					 dfr=dfr , BIFIEobj=BIFIEobj )				
	dfr <- clean_summary_table( dfr=dfr , RR=RR , se=se , Nimp=Nimp )	
	
				
	# create vector of parameter names
	parnames <- base::paste0( dfr$parm   , "_" , dfr$groupvar , dfr$groupval )


	#@@@@***
	# multiple groupings
	dfr <- BIFIE_table_multiple_groupings( dfr , res00 )
	#@@@@***
						
	
	#*************************** OUTPUT ***************************************
	s2 <- base::Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- base::list( "stat" = dfr , 
			"output" = res , 	"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac , "GG"=GG ,			
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"parnames" = parnames , "CALL"= cl)
	base::class(res1) <- "BIFIE.by"
	base::return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.by function
summary.BIFIE.by <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	base::cat("Statistical Inference for User Defined Function \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}
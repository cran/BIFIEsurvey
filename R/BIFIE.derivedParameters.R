
####################################################
# statistical inference for derived parameters
BIFIE.derivedParameters <- function( BIFIE.method , derived.parameters , type=NULL ){
	
	cl <- match.call()	
	s1 <- Sys.time()		
	
	object <- res1 <- BIFIE.method		
	parnames <- res1$parnames
			
	# extract replicated parameters
	parsres <- extract.replicated.pars( BIFIE.method = res1 , type=type )

	#*****			
	pars0 <- parsres$parsM
	pars0.rep <- parsres$parsrepM	
	rownames(pars0) <- parnames		
	rownames(pars0.rep) <- parnames			
	
	###########################
	allformulas <- derived.parameters[[1]]
	FF <- length(derived.parameters)
	if (FF>1){
	for (ff in 2:FF){
		t1 <- terms( allformulas)    
		t2 <- paste( c( attr( t1 , "term.labels" ) , 
			attr(  terms( derived.parameters[[ff]] ) , "term.labels" )  ),
			collapse= " + " )
		allformulas  <- as.formula( paste( " ~ 0 + " , t2 ) )
				}
				}
	# create matrices of derived parameters
	der.pars <- model.matrix( allformulas , as.data.frame( t(pars0) ) )
	colnames(der.pars) <- names(derived.parameters)
	der.pars.rep <- model.matrix( allformulas , as.data.frame( t(pars0.rep) ) )
	colnames(der.pars.rep) <- names(derived.parameters)
	fayfac <- res1$fayfac
	NP <- ncol(der.pars)
	Cdes <- matrix( 1 , ncol=NP , nrow=1 )
	Ccols <- which( colSums( abs( Cdes) ) > 0 )
	parsM <- as.matrix( t( der.pars ) )
	parsrepM <- as.matrix( t( der.pars.rep ) )
	rdes <- c(0)
	# compute covariance matrices
	res0 <- .Call("bifie_waldtest" ,  parsM = parsM , parsrepM = parsrepM , 
						Cdes , rdes , Ccols - 1 , fayfac=fayfac ,
						PACKAGE="BIFIEsurvey")
	var_w <- res0$var_w
	var_b <- res0$var_b
	Nimp <- res1$Nimp
	# total variance
	var_tot <- var_w  + ( 1 + 1/Nimp ) * var_b 
	parmlabel <- names(derived.parameters)
	# parameters and standard errors
	stat <- data.frame( "parmlabel" = parmlabel , "coef"= rowMeans( parsM ) ,
				"se" = sqrt( diag( var_tot ) ) ) 
	# pars_fmi[pp] = ( 1.0 + 1/Nimp2) * pars_varBetween[pp] / pow(pars_se[pp] + eps,2.0) ;
	eps <- 1E-10
	
	stat$t <- stat$coef / stat$se
	stat$p <- 2* pnorm( - abs( stat$t ) )
	stat$fmi <-  ( 1+1/Nimp) * diag(var_b) / ( stat$se^2 + eps )           
	stat$VarMI <- diag( var_b )
	stat$VarRep <- diag( var_w )               
	rownames(stat) <- NULL         
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )	
	res <- list( "stat" = stat , "coef"= rowMeans( parsM ) ,
				"se" = sqrt( diag( var_tot ) ) ,
				"vcov" = var_tot , "Nimp" = Nimp , "fayfac" = fayfac ,
				"N"=res1$N , "RR"=res1$RR , 
				"allformulas" = allformulas , "CALL"=cl ,
				"timediff" = timediff 	,
				"derived.parameters" = derived.parameters ,
				"parnames" = names(derived.parameters)
						)	
	class(res) <- "BIFIE.derivedParameters"
	return(res)
	}
#########################################################

####################################################################################
# summary for BIFIE.derivedParameters function
summary.BIFIE.derivedParameters <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Formulas for Derived Parameters \n")	
	FF <- length( object$derived.parameters)
	for (ff in 1:FF){
		# ff <- 1
		# cat( paste0( parnames[ff]  , paste( derived.parameters[[ff]] ) ) , "\n") 
		cat( paste0(  object$parnames[ff] , " := " , 
					attr(  terms( object$derived.parameters[[ff]] ) , "term.labels" )  , 
						collapse=" " ) , 
								"\n") 	
						}
	#***
	cat("\nStatistical Inference for Derived Parameters \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}	
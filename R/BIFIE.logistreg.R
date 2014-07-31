

#######################################################################
# Linear regression
BIFIE.logistreg <- function( BIFIEobj , dep=NULL , pre=NULL  , 
            formula=NULL ,   
			group=NULL , group_values=NULL , se=TRUE , eps = 1E-8 , maxiter=100){
	#****
	s1 <- Sys.time()
	bifieobj <- BIFIEobj
	if (bifieobj$cdata){
	    fomula_vars <- NULL
		if (! is.null(formula) ){
			formula_vars <- all.vars( formula )
							}
		varnames <- unique( c( dep , pre , group , "one" , formula_vars ) )
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

	#*** look for formula objects
	if ( ! is.null( formula) ){
	    cat("|*** Data Preparation ") ; flush.console()
	    bifieobj2 <- datalistM
		colnames(bifieobj2) <- varnames
		if ( is.null(group) ){ group <- "one" ; group_values <- 1 }
		bifieobj2 <- as.data.frame( bifieobj2 )
	    m1 <- model.matrix(formula , data=bifieobj2)
		#***
        m0 <- m1
		m1 <- matrix( NA , nrow=nrow(bifieobj2) , ncol=ncol(m0) )
		m1[ match( rownames(m0),rownames(bifieobj2) ) , ] <- m0
		colnames(m1) <- colnames(m0)
		#****
		dep <- rownames( attr( terms(formula) ,"factors") )[1]
		pre <- colnames( m1 )
		datalistM <- as.matrix( cbind( bifieobj2[ , dep  ] , m1 , bifieobj2[,group] ) )
		varnames <- c( dep , pre , group )	
        cat("\n")		
		}
	if (RR==1){ RR <- 0 }
	if ( ! se ){ 
		wgtrep <- matrix( wgt , ncol=1 )
		RR <- 0
				}	
	
	dep_index <- unlist( sapply( dep , FUN = function(vv){ 
                        which( varnames == vv ) } ) )
	pre_index <- unlist( sapply( pre , FUN = function(vv){ 
                        which( varnames == vv ) } ) )
	# vars values
	VV <- length(pre)
					
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
	res <- .Call("bifie_logistreg" ,  datalistM , wgt_ , as.matrix(wgtrep) , dep_index -1 , 
				pre_index - 1 ,  fayfac ,    Nimp ,  group_index -  1, group_values ,
				eps , maxiter , PACKAGE="BIFIEsurvey" )

		
	GG <- length(group_values)
#	ZZ <- nrow(itempair_index )	
    ZZ <- VV
	p1 <- c( rep("b",VV) )
    p2 <- c( pre )
	dfr <- data.frame( "parameter" = rep(p1,GG)	)
	dfr$var <- rep(p2,GG)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , each=ZZ )
	             }							 
	dfr$Ncases <- rep( rowMeans( res$ncasesM ) , each=ZZ )
	dfr$Nweight <- rep( rowMeans( res$sumwgtM ) , each=ZZ )	
	
		
	dfr$est <- res$regrcoefL$pars
	dfr$SE <- res$regrcoefL$pars_se
	dfr$t <- round( dfr$est / dfr$SE , 2 )
	dfr$p <- pnorm( - abs( dfr$t ) ) * 2
	dfr$fmi <- res$regrcoefL$pars_fmi
	dfr$VarMI <- res$regrcoefL$pars_varBetween
	dfr$VarRep <- res$regrcoefL$pars_varWithin
	if ( ( ! se ) &  ( RR==0 ) ){					
		dfr$t <- dfr$p <- dfr$SE <- dfr$fmi <- dfr$VarMI <- dfr$VarRep <- NULL
				}				

	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )
	parnames <- paste0( dfr$parameter   , "_" , dfr$var , 
			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
			ifelse( ! nogroupL , dfr$groupval , "" ) )
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr ,
			"output" = res , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"GG"=GG , "parnames" = parnames)
	class(res1) <- "BIFIE.logistreg"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.linreg function
summary.BIFIE.logistreg <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Statistical Inference for Logistic Regression \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}


#######################################################################
# BIFIE univar.test
BIFIE.univar.test <- function( BIFIE.method ){
	#****
	s1 <- Sys.time()	
	res5 <- BIFIE.method	
	if (res5$group == "one"){ stop("This function can only be applied with a grouping variable.\n")}
		
	mean1M <- res5$output$mean1M
	sd1M <- res5$output$sd1M
	sumweightM <- res5$output$sumweightM
	GG <- res5$GG
	group_values <- ( res5$stat$groupval )[1:GG]
	mean1repM <- res5$output$mean1repM
	sd1repM <- res5$output$sd1repM
	sumweightrepM <- res5$output$sumweightrepM
	fayfac <- res5$fayfac
	VV <- res5$VV
	vars <- res5$vars
    RR <- res5$RR
	if (RR==1){ RR <- 0 }
	group <- res5$group
	N <- res5$N
	Nimp <- res5$Nimp
		
	#*****
	# Rcpp call
	res <- .Call("bifie_test_univar" , mean1M , sd1M , sumweightM , GG , group_values ,
			mean1repM , sd1repM , sumweightrepM  , fayfac , PACKAGE="BIFIEsurvey")

#	res <- bifie_test_univar(  mean1M , sd1M , sumweightM , GG , group_values ,
#			mean1repM , sd1repM , sumweightrepM  , fayfac )

			
	#****
	# output eta^2		
	dfr <- data.frame( "var" = vars	, "group" = group  )	
    dfr$eta2 <- res$eta2L$pars^2	
	dfr$eta <- res$eta2L$pars
	dfr$eta_SE <- res$eta2L$pars_se
	dfr$fmi <- res$eta2L$pars_fmi
	dfr$VarMI <- res$eta2L$pars_varBetween
	dfr$VarRep <- res$eta2L$pars_varWithin
	if (RR==0){				
		dfr$SE <- dfr$fmi <- dfr$VarMI <- dfr$VarRep <- NULL
				}				
	stat.eta2 <- dfr
	
	#****
	# output d values
	group_values_matrix <- res$group_values_matrix	
	ZZ <- nrow(group_values_matrix)
	dfr <- data.frame( "var" = rep(vars,each=ZZ)	)	
	g2 <- group_values_matrix[ rep(1:ZZ , VV) , , drop=FALSE]
	g2 <- data.frame(g2)
	colnames(g2) <- c("groupval1" , "groupval2")
	dfr$group <- group
	dfr <- cbind( dfr , g2 )
	r5 <- res5$stat
	h1 <- NULL
	h4 <- h3 <- h2 <- NULL
    for (kk in 1:2){
		group_values_matrix[,kk] <- match( group_values_matrix[,kk] , group_values )
		}
	for (vv in 1:VV){
		h1 <- c(h1 , r5[ (vv-1)*GG + group_values_matrix[,1] , "M" ] )
		h2 <- c(h2 , r5[ (vv-1)*GG + group_values_matrix[,2] , "M" ] )		
		h3 <- c(h3, r5[ (vv-1)*GG + group_values_matrix[,1] , "SD" ] )
		h4 <- c(h4 , r5[ (vv-1)*GG + group_values_matrix[,2] , "SD" ] )		
					}
	dfr$M1 <- h1
    dfr$M2 <- h2
    dfr$SD <- sqrt( ( h3^2 + h4^2 ) / 2 )
	dfr$d <- res$dstatL$pars
	dfr$d_SE <- res$dstatL$pars_se
	dfr$t <- dfr$d / dfr$d_SE
	dfr$p <- pnorm(  - abs(dfr$t ) )*2
	dfr$fmi <- res$dstatL$pars_fmi
	dfr$VarMI <- res$dstatL$pars_varBetween
	dfr$VarRep <- res$dstatL$pars_varWithin
	if (RR==0){				
		dfr$SE <- dfr$fmi <- dfr$VarMI <- dfr$VarRep <- NULL
				}				
	stat.dstat <- dfr	
	
	#*****
	# F statistics

	dfr <- NULL
	for (vv in 1:VV){
		Cdes <- matrix( 0 , nrow=GG-1 , ncol=GG*VV )
		indvec <- 1:(GG-1)
		for (zz in indvec ){
			Cdes[ zz , c(zz + (vv-1)*GG ,zz+1 + (vv-1)*GG ) ] <- c(1,-1)
					}
		rdes <- rep(0,GG-1)
		wres5 <- BIFIE.waldtest( res5 , Cdes=Cdes , rdes=rdes )		
		dfr <- rbind( dfr , wres5$stat.D )
				}
	dfr <- data.frame( "variable" = vars , "group" = group , dfr )	
	stat.F <- dfr
	
	
	
	parnames <- NULL
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.F"= stat.F , "stat.eta" = stat.eta2 , "stat.dstat" = stat.dstat , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"GG"=GG , "parnames" = parnames)
	class(res1) <- "BIFIE.univar.test"
	return(res1)
		}
###################################################################################

####################################################################################
# summary
summary.BIFIE.univar.test <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("F Test (ANOVA) \n")	
	obji <- object$stat.F
	print.object.summary( obji , digits=digits )	
	cat("\nEta Squared \n")	
	obji <- object$stat.eta
	print.object.summary( obji , digits=digits )
	cat("\nCohen's d Statistic \n")	
	obji <- object$stat.dstat
	print.object.summary( obji , digits=digits )	
			}
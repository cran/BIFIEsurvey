###########################################################
# BIFIE.data objects for designs with jackknife zones
BIFIE.data.jack <- function( data , wgt=NULL , pv_vars = NULL ,
	jkzone=NULL , jkrep=NULL , jktype="JK_TIMSS" ,
	jkfac=NULL , fayfac = 1
		){
	
    data <- as.data.frame( data )
	
	#**** defaults for TIMSS	
	if (jktype == "JK_TIMSS"){
          jkrep <- "JKREP"
          jkzone <- "JKZONE"
		  wgt <- "TOTWGT"
          jkfac <- 2
				}
				
	#******** generate replicate weights
	RR <- max( data[,jkzone] )
#	datarep <- matrix( NA , nrow=nrow(data) , ncol=RR )
#	colnames(datarep) <- paste0("R",1:RR)
	prblen <- 10
	prbar <- BIFIE.progressbar( ops = RR , prblen = prblen )
	cat("+++ Generate replicate weights\n")
	cat(paste0("|" , paste0(rep("*",prblen), collapse="") , "|\n|")) ; flush.console()
	
	addname <- 10^( floor( log( RR+.5 , 10 ) )  + 1 )
#	colnames(datarep) <-  paste0("w_fstr" , substring( paste0(addname +1:RR),2) )	
#	for (rr in 1:RR){    # rr <- 1
#		vv <- paste0("w_fstr" , substring( addname +rr,2) )
#		datarep[ , rr ] <- round( ifelse( as.numeric(paste0(data[,jkzone])) == rr , 
#							jkfac*data[,wgt] * data[,jkrep] , 1*data[,wgt] ) , 5 )
#		colnames(datarep)[rr] <- vv
#		if ( prbar[rr] == 1 ){ cat("-" ) ; flush.console() }
#				}

	datarep <- .Call( "bifie_jack_timss" , wgt_= data[,wgt] , data[,jkzone]-1 , data[,jkrep] , RR , jkfac ,
					prbar , PACKAGE="BIFIEsurvey" )

	colnames(datarep) <- paste0("w_fstr" , substring( paste0(addname +1:RR),2) )	
	
	cat("|\n")
	#******** generate replicated datasets for datasets
	if ( is.null( pv_vars) ){ datalist <- data  }
	if ( ! is.null( pv_vars )){
		dfr <- NULL
		VV <- length(pv_vars)
		for (vv in 1:VV){
			vv1 <- pv_vars[vv]
			ind.vv1 <- which( substring( colnames(data) , 1 , nchar( vv1 ) ) == pv_vars[vv] )
			Nimp <- length(ind.vv1)
			dfr2 <- data.frame( "variable" = vv1 , "var_index" = vv , "data_index" = ind.vv1 ,
							 "impdata_index"=1:Nimp ) 			
			dfr <- rbind( dfr , dfr2 )
			}
		sel_ind <- setdiff( 1:( ncol(data) ) , dfr$data_index )
		data0 <- data[ , sel_ind ]
		V0 <- ncol(data0)
		newvars <- seq( V0+1 , V0+VV )
		datalist <- as.list( 1:Nimp )
		for (ii in 1:Nimp ){
			dat1 <- data.frame( data0 , data[ , dfr[ dfr$impdata_index == ii  , "data_index" ] ] )
			colnames(dat1)[ newvars ] <- pv_vars
			datalist[[ii]] <- dat1 
						}
					}
	#*** create BIFIE.data object
	bifiedat <- BIFIE.data( datalist , wgt = data[, wgt ] , wgtrep = datarep , fayfac = fayfac )
	return(bifiedat)
}
###############################################################################
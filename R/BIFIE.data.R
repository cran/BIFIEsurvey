
##################################################################
# convert a list of multiply imputed datasets into an object
# of class BIFIEdata
BIFIE.data <- function( data.list , wgt=NULL , wgtrep=NULL , fayfac=1){
	if ( ( is.list( data.list ) ) & ( is.data.frame( data.list) ) ){ 
	    h1 <- data.list
		data.list <- list( 1 )
		data.list[[1]] <- h1
				}
    FF <- length( data.list)
    Nimp <- FF
    N <- nrow( data.list[[1]] )
    V <- ncol( data.list[[1]] )
    dat1 <- data.list[[1]]
	cn <- c( colnames(dat1) , "one" )
    N <- nrow(dat1)    
    p1 <- sapply( 1:V , FUN = function(vv){ is.numeric( dat1[,vv] ) } )
    notnum <- which( ! p1 )
    datalistM <- matrix( NA , nrow=N*Nimp , V + 1)
	cat("+++ Generate BIFIE.data object\n")
	cat(paste0( "|" , paste0( rep("*" , FF) , collapse="") , "|\n|" ))
	#****
	# weights
	if ( is.null(wgt) ){ wgt <- rep(1,N) }	
	wgt <- as.numeric( wgt )
	if ( is.null(wgtrep) ){ wgtrep <- matrix( wgt , nrow=N , ncol=1 ) }
	wgtrep <- as.matrix( wgtrep )
    for (ff in 1:FF){
        dat1 <- data.list[[ff]] 
        for (vv in notnum){ 
			dat1[,vv] <- as.numeric(paste( dat1[,notnum] ) ) 			
		}
		dat1$one <- 1		
        dat1 <- as.matrix( dat1) 
        datalistM[ 1:N + N*(ff-1) , ] <- dat1 
		cat("-") ; flush.console()
                }
	cat("|\n")	
    res <- list( "datalistM" = datalistM , "wgt" = wgt , "wgtrep" = wgtrep ,
        "Nimp" = Nimp , "N"= N , "dat1" = dat1  , "varnames" = cn , "fayfac"= fayfac ,
		"RR"= ncol(wgtrep) , "time" = rep( Sys.time(),2) )
    class(res) <- "BIFIEdata"
    return(res)
    }            
########################################################################
summary.BIFIEdata <- function(object , ... ){
	cat("------------------------------------------------------------\n")
	d1 <- packageDescription("BIFIEsurvey")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( paste0("Function '" , class(object) ) )
	cat("'\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
    cat( "Number of persons =" , object$N , "\n" )    
    cat( "Number of imputed datasets =" , object$Nimp , "\n" ) 
    cat( "Number of Jackknife zones per dataset =" , object$RR , "\n" ) 
    cat( "Fay factor =" , round( object$fayfac , 5 ) , "\n\n" ) 	
						}
######################################################################						
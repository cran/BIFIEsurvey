
##################################################################
# Convert a list of multiply imputed datasets into an object
#      of class BIFIEdata
BIFIE.data <- function( data.list , wgt=NULL , wgtrep=NULL , fayfac=1 , cdata=FALSE){
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
	if ( is.character(wgt) & ( length(wgt) == 1 ) ){
		wgt <- data.list[[1]][ , wgt ]
				}
	
	if ( is.null(wgt) ){ wgt <- rep(1,N) }	
	wgt <- as.numeric( wgt )
	if ( is.null(wgtrep) ){ wgtrep <- matrix( wgt , nrow=N , ncol=1 ) }
	wgtrep <- as.matrix( wgtrep )
    for (ff in 1:FF){  # imputed dataset ff
        dat1 <- data.list[[ff]] 
        for (vv in notnum){ 		
#			dat1[,vv] <- as.numeric(paste( dat1[,vv] ) ) 			
			dat1[,vv] <- as.numeric( dat1[,vv] ) 			
		}
		dat1$one <- 1		
        dat1 <- as.matrix( dat1) 	
        datalistM[ 1:N + N*(ff-1) , ] <- dat1 
		cat("-") ; flush.console()
                }
	cat("|\n")	
    res <- list( "datalistM" = datalistM , "wgt" = wgt , "wgtrep" = wgtrep ,
        "Nimp" = Nimp , "N"= N , "dat1" = dat1  , "varnames" = cn , "fayfac"= fayfac ,
		"RR"= ncol(wgtrep) , "time" = Sys.time() )
	res$cdata <- FALSE	
    class(res) <- "BIFIEdata"	
	#***** variable names and transformations
	VV <- length(res$varnames)
	res$Nvars <- VV
	dfr2 <- data.frame( "index" = 1:VV , "variable" = res$varnames ,
				"variable_orig" = res$varnames , "source"="indata")
	res$variables <- dfr2
	if ( cdata ){
		res <- BIFIE.BIFIEdata2BIFIEcdata( bifieobj =res , varnames=NULL )		
				}
    return(res)
    }            
########################################################################
summary.BIFIEdata <- function(object , ... ){
	cat("------------------------------------------------------------\n")
	d1 <- packageDescription("BIFIEsurvey")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( paste0("Function '" , class(object) ) )
	cat("'" )		
	if( object$cdata ){ cat("\nCompact BIFIEdata") }
	cat("\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[1] ) , "\n" )
    cat( "Number of persons =" , object$N , "\n" )    
    cat( "Number of variables =" , object$Nvars , "\n" ) 	
    cat( "Number of imputed datasets =" , object$Nimp , "\n" ) 
    cat( "Number of Jackknife zones per dataset =" , object$RR , "\n" ) 
    cat( "Fay factor =" , round( object$fayfac , 5 ) , "\n\n" ) 	
	x2 <- BIFIE_object_size(object)
	cat( "Object size: \n  " ) 			
	cat( "Total: " , x2$value_string , "\n")
	
	obj1 <- "datalistM"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
    
	obj1 <- "datalistM_ind"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")

	obj1 <- "datalistM_imputed"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
	
	obj1 <- "datalistM_impindex"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")
	
	obj1 <- "dat1"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")	

	obj1 <- "wgtrep"
	x2 <- BIFIE_object_size(object[[ obj1 ]] )
	cat( paste0( "   $" , obj1 , " :" ) , x2$value_string , "\n")	

# Revalpr("x2")	
#    x1 <- utils::object.size(object)
#    cat( "Object size: " , round( as.numeric(x1 / 1024^2) , 3 ) , "MB" ) 		
#    cat( " | " , round( as.numeric(x1 / 1024^3) , 5 ) , "GB\n" ) 			
						}
######################################################################	

##################################################
# output object sizes in MB and GB
BIFIE_object_size <- function( x1 ){
	x1 <- utils::object.size(x1)
	vals <- c( x1 , as.numeric(x1 / 1024^2) ,
					as.numeric(x1 / 1024^3)  )
	names(vals) <- c("B" , "MB" , "GB" )		
    res <- list( "value" = vals )
	res$value_string <- paste0( round( vals[2] , 3 ) , 
			" MB | " , round( vals[3] , 5 ) ,	" GB" ) 		
    return(res)
			}
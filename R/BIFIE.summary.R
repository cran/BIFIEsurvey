
BIFIE.summary <- function(object , print.time=TRUE){
	cat("------------------------------------------------------------\n")
	d1 <- packageDescription("BIFIEsurvey")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
	cat( paste0("Function '" , class(object) ) )
	if ( class(object) == "BIFIE.waldtest" ){
	    cat( paste0( "' for BIFIE method '" , object$class.BIFIE.method ) )
					}							
	cat("'" )	
	cat("\n\nCall:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
	# cat("\n\n")
	cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
	if (print.time){
		cat("Computation time:" , print(object$time[2] - object$time[1] ), "\n\n") 
			} else { 
		cat("\n")
			}
    cat( "Number of persons =" , object$N , "\n" )    
    cat( "Number of imputed datasets =" , object$Nimp , "\n" ) 
    cat( "Number of Jackknife zones per dataset =" , object$RR , "\n" ) 
    cat( "Fay factor =" , round( object$fayfac , 5 ) , "\n\n" ) 	
						}
## File Name: BIFIE.data.R
## File Version: 1.473


# Convert a list of multiply imputed datasets into an object of class BIFIEdata
BIFIE.data <- function( data.list, wgt=NULL, wgtrep=NULL, fayfac=1,
        pv_vars=NULL, pvpre=NULL, cdata=FALSE, NMI=FALSE )
{

    cl <- match.call()

    #**** handling pv_vars
    if ( ! is.null(pv_vars) ){
        if ( is.null(pvpre) ){
            jktype <- "JK_TIMSS"
        } else {
            jktype <- "RW_PISA"
        }
        if (!is.null(pvpre)){
            cn_data <- colnames(data.list)
            pv_vars <- BIFIE_data_select_pv_vars(pvpre=pvpre, cn_data=cn_data)
        }
        data.list <- BIFIE_data_pv_vars_create_datlist(pvpre=pvpre, pv_vars=pv_vars,
                jktype=jktype, data=data.list)
    }

    # subroutine for preparation of nested multiple imputations
    res0 <- BIFIE_data_nested_MI( data.list=data.list, NMI=NMI )
    data.list <- res0$data.list
    Nimp_NMI <- res0$Nimp_NMI

    if ( ( is.list( data.list ) ) & ( is.data.frame( data.list) ) ){
        h1 <- data.list
        data.list <- list( 1 )
        data.list[[1]] <- h1
    }

    FF <- length( data.list)
    Nimp <- FF
    if ( sum( colnames(data.list[[1]]) %in% "one" ) > 0 ){
        cat("Variable 'one' in datasets is replaced by a constant variable")
        cat(" containing only ones!\n" )
        for (ii in 1:Nimp){
            data.list[[ii]][, "one"] <- NULL
        }
    }
    N <- nrow( data.list[[1]] )
    V <- ncol( data.list[[1]] )
    dat1 <- data.list[[1]]

    cn <- c( colnames(dat1), "one" )
    N <- nrow(dat1)

    p1 <- sapply( 1:V, FUN=function(vv){ is.numeric( dat1[,vv] ) } )
    notnum <- which( ! p1 )
    datalistM <- matrix( NA, nrow=N*Nimp, V + 1)
    cat("+++ Generate BIFIE.data object\n")
    cat(paste0( "|", paste0( rep("*", FF), collapse=""), "|\n|" ))
    #****
    # weights
    if ( is.character(wgt) & ( length(wgt)==1 ) ){
        wgt <- data.list[[1]][, wgt ]
    }

    if ( is.null(wgt) ){ wgt <- rep(1,N) }
    wgt <- as.numeric( wgt )
    if ( is.null(wgtrep) ){ wgtrep <- matrix( wgt, nrow=N, ncol=1 ) }
    wgtrep <- as.matrix( wgtrep )
    for (ff in 1:FF){  # imputed dataset ff
        dat1 <- data.list[[ff]]
        for (vv in notnum){
            dat1[,vv] <- as.numeric( dat1[,vv] )
        }
        dat1$one <- 1
        dat1 <- as.matrix( dat1)
        datalistM[ 1:N + N*(ff-1), ] <- dat1
        cat("-") ; flush.console()
    }
    cat("|\n")
    wgtrep <- as.matrix(wgtrep)
    res <- list( "datalistM"=datalistM, "wgt"=wgt, "wgtrep"=wgtrep,
        "Nimp"=Nimp, "N"=N, "dat1"=dat1, "varnames"=cn, "fayfac"=fayfac,
        "RR"=ncol(wgtrep), "time"=Sys.time(), "CALL"=cl )
    res$NMI <- NMI
    res$Nimp_NMI <- Nimp_NMI
    res$cdata <- FALSE
    class(res) <- "BIFIEdata"
    #***** variable names and transformations
    VV <- length(res$varnames)
    res$Nvars <- VV
    dfr2 <- data.frame( "index"=1:VV, "variable"=res$varnames,
                "variable_orig"=res$varnames, "source"="indata")
    res$variables <- dfr2
    if ( cdata ){
        res <- BIFIE.BIFIEdata2BIFIEcdata( bifieobj=res, varnames=NULL )
                }
    return(res)
}

#**************** print method ***********************
print.BIFIEdata <- function(x,...){
    cat("Object of class 'BIFIEdata'\nCall: ")
    print( x$CALL )
    #*** multiply imputed data
    if ( ! x$NMI ){
        cat("MI data with", x$Nimp,"datasets\n")
    }
    #*** nested multiply imputed data
    if ( x$NMI ){
        v1 <- paste0( "NMI data with ", x$Nimp_NMI[1]," between datasets and ",
              x$Nimp_NMI[2], " within datasets\n")
        cat(v1)
    }
    v1 <- paste0( x$RR, " replication weights with fayfac=",
                    round(x$fayfac,3), " \n" )
    cat(v1)
    v1 <- paste0( x$N, " cases and ",    x$Nvars, " variables \n" )
    cat(v1)
}

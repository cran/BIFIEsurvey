\name{BIFIE.by}
\alias{BIFIE.by}
\alias{summary.BIFIE.by}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Statistics for User Defined Functions
}

\description{
Computes statistics for user defined functions.
}

\usage{
BIFIE.by( BIFIEobj , vars , userfct , userparnames=NULL ,
     group=NULL , group_values=NULL , se=TRUE )

\method{summary}{BIFIE.by}(object,digits=4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars}{
Vector of variables for which statistics should be computed
}
\item{userfct}{User defined function. This function must include
a matrix \code{X} and a weight vector \code{w} as
arguments. The value of this function must be a vector.
}
\item{userparnames}{An optional vector of parameter names for the
	value of \code{userfct}.	
}
  \item{group}{
Optional grouping variable
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{object}{Object of class \code{BIFIE.by}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with statistics defined in \code{userfct}}
 \item{output}{Extensive output with all replicated statistics}
 \item{\dots}{More values}
}

%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survey:svyby]{svyby}} (\pkg{survey})
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bifieobj <- BIFIE.data( data.list=data.timss1 , wgt=  data.timss1[[1]]$TOTWGT ,
           wgtrep=data.timssrep[, -1 ] )

#*** Model 1: Weighted means (as a toy example)
userfct <- function(X,w){
        pars <- c( weighted.mean( X[,1] , w ) , 
                 weighted.mean(X[,2] , w )   )
        return(pars)
                        }
res1 <-  BIFIE.by( bifieobj , vars =c("ASMMAT" , "migrant" , "books" ) , 
                userfct=userfct , userparnames =c("MW_MAT", "MW_Migr") , 
                group="female" )
summary(res1)

\dontrun{
#**** Model 2: Robust linear model

# (1) start from scratch to formulate the user function for X and w
dat1 <- bifieobj$dat1
vars <- c("ASMMAT" , "migrant" , "books" )
X <- dat1[,vars]
w <- bifieobj$wgt
library(MASS)
# ASMMAT ~ migrant + books
mod <- MASS::rlm( X[,1] ~  as.matrix( X[, -1 ] ) , weights = w )
coef(mod)
# (2) define a user function "my_rlm"
my_rlm <- function(X,w){
    mod <- MASS::rlm( X[,1] ~  as.matrix( X[, -1 ] ) , weights = w )
    return( coef(mod) )
                }
# (3) estimate model
res2 <-  BIFIE.by( bifieobj , vars  , userfct=my_rlm , group="female" , group_values=0:1)
summary(res2)
# estimate model without computing standard errors
res2a <-  BIFIE.by( bifieobj , vars  , userfct=my_rlm , group="female" , se=FALSE)
summary(res2a)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{User defined function}
\keyword{summary}

\name{BIFIE.data}
\alias{BIFIE.data}
\alias{summary.BIFIEdata}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Creates an Object of Class \code{BIFIEdata}
}

\description{
This function creates an object of class \code{BIFIEdata}.
}
\usage{
BIFIE.data(data.list, wgt = NULL, wgtrep = NULL , fayfac = 1, cdata=FALSE)

\method{summary}{BIFIEdata}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.list}{
List of multiply imputed datasets
}
  \item{wgt}{
A string indicating the label of case weight or a vector
containing all case weights.
}
  \item{wgtrep}{
Optional vector of replicate weights
}
\item{fayfac}{Fay factor for calculating standard errors}
\item{cdata}{An optional logical indicating whether the \code{BIFIEdata}
object should be compactly saved. The default is \code{FALSE}.}
\item{object}{Object of class \code{BIFIEdata}}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An object of class \code{BIFIEdata} saved in a non-compact
or compact way, see value \code{cdata}. The following entries are
included in the list:

\item{datalistM}{Stacked list of imputed datasets (if \code{cdata=FALSE})}
\item{wgt}{Vector with case weights}
\item{wgtrep}{Matrix with replicate weights}
\item{Nimp}{Number of imputed datasets}
\item{N}{Number of observations in a dataset}
\item{dat1}{Last imputed dataset}
\item{varnames}{Vector with variable names}
\item{fayfac}{Fay factor}
\item{RR}{Number of replicate weights}
\item{cdata}{Logical indicating whether the \code{BIFIEdata} object is in
	compact format (\code{cdata=TRUE}) or in a non-compact format
	(\code{cdata=FALSE}).
}
\item{Nvars}{Number of variables}
\item{variables}{Data frame including some informations about variables.
All transformations are saved in the column \code{source}.}
\item{datalistM_ind}{Data frame with response indicators 
(if \code{cdata=TRUE})}
\item{datalistM_imputed}{Data frame with imputed values 
(if \code{cdata=TRUE})}		
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
See \code{\link{BIFIE.data.transform}} for data transformations on
\code{BIFIEdata} objects.

For saving and loading \code{BIFIEdata} objects see 
\code{\link{save.BIFIEdata}}.

For converting PIRLS/TIMSS or PISA datasets into \code{BIFIEdata}
objects see \code{\link{BIFIE.data.jack}}.
}

\examples{
#############################################################################
# EXAMPLE 1: Create BIFIEdata object with multiply imputed TIMSS data
#############################################################################
data(data.timss1)
data(data.timssrep)

bdat <- BIFIE.data( data.list=data.timss1 , wgt= data.timss1[[1]]$TOTWGT ,
            wgtrep=data.timssrep[, -1 ] )
summary(bdat)            
# create BIFIEdata object in a compact way
bdat2 <- BIFIE.data( data.list=data.timss1 , wgt= data.timss1[[1]]$TOTWGT ,
            wgtrep=data.timssrep[, -1 ] , cdata=TRUE)
summary(bdat2)  

\dontrun{
#############################################################################
# EXAMPLE 2: Create BIFIEdata object with one dataset
#############################################################################
data(data.timss2)

# use first dataset with missing data from data.timss2
bdat <- BIFIE.data( data.list=data.timss2[[1]] , wgt=data.timss2[[1]]$TOTWGT )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
\keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

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
BIFIE.data(data.list, wgt = NULL, wgtrep = NULL , fayfac = 1)

\method{summary}{BIFIEdata}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.list}{
List of multiply imputed datasets
}
  \item{wgt}{
Optional vector of person weights
}
  \item{wgtrep}{
Optional vector of replicate weights
}
\item{fayfac}{Fay factor for calculating standard errors}
\item{object}{Object of class \code{BIFIEdata}}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Object of class \code{BIFIEdata}
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

% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
% }

\examples{
#############################################################################
# EXAMPLE 1: Create BIFIEdata object with multiply imputed TIMSS data
#############################################################################
data(data.timss1)
data(data.timssrep)

bdat <- BIFIE.data( data.list=data.timss1 , wgt= data.timss1[[1]]$TOTWGT ,
            wgtrep=data.timssrep[, -1 ] )
summary(bdat)            

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

\name{BIFIE.freq}
\alias{BIFIE.freq}
\alias{summary.BIFIE.freq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Frequency Statistics
}

\description{
Computes some absolute and relative frequencies.
}

\usage{
BIFIE.freq(BIFIEobj, vars, group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.freq}(object,digits=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars}{
Vector of variables for which statistics should be computed
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
\item{object}{Object of class \code{BIFIE.univar}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with frequency statistics}
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
\code{\link[survey:svytable]{svytable}} (\pkg{survey}),
\code{\link[intsvy:timss.table]{timss.table}} (\pkg{intsvy}),
\code{\link[Hmisc:wtd.stats]{wtd.table}} (\pkg{Hmisc}),
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIE.data( data.list=data.timss1 , wgt=  data.timss1[[1]]$TOTWGT ,
           wgtrep=data.timssrep[, -1 ] )

# Frequencies for three variables
res1 <- BIFIE.freq( bdat , vars=c("lang" , "books" , "migrant" )  )
summary(res1)

# Frequencies splitted by gender
res2 <- BIFIE.freq( bdat , vars=c("lang" , "books" , "migrant" ) , 
              group="female" , group_values=0:1 )
summary(res2)                    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}
\name{BIFIE.correl}
\alias{BIFIE.correl}
\alias{summary.BIFIE.correl}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Correlations and Covariances
}

\description{
Computes correlations and covariances
}

\usage{
BIFIE.correl(BIFIEobj, vars, group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.correl}(object,digits=4,...)
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
 \item{stat.cor}{Data frame with correlation statistics}
 \item{stat.cov}{Data frame with covariance statistics}
 \item{cor_matrix}{List of estimated correlation matrices}
 \item{cov_matrix}{List of estimated covariance matrices} 
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
\code{\link[stats:cov.wt]{cov.wt}} (\pkg{stats}),	
\code{\link[intsvy:timss.rho]{timss.rho}} (\pkg{intsvy}),
\code{\link[intsvy:timss.rho.pv]{timss.rho.pv}} (\pkg{intsvy}),
\code{\link[Hmisc:rcorr]{rcorr}} (\pkg{Hmisc}),
\code{\link[miceadds:ma.wtd.statNA]{ma.wtd.corNA}} (\pkg{miceadds})
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

# Correlations splitted by gender
res1 <- BIFIE.correl( bdat , vars=c("lang" , "books" , "migrant" ) , 
              group="female" , group_values=0:1 )
summary(res1)  

# Correlations splitted by gender: no statistical inference (se=FALSE)
res1a <- BIFIE.correl( bdat , vars=c("lang" , "books" , "migrant" ) , 
              group="female" , group_values=0:1 , se=FALSE)
summary(res1a)  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}

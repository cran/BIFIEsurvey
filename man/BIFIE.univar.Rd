\name{BIFIE.univar}
\alias{BIFIE.univar}
\alias{summary.BIFIE.univar}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Univariate Descriptive Statistics (Means and Standard Deviations)
}

\description{
Computes some univariate descriptive statistics 
(means and standard deviations).
}

\usage{
BIFIE.univar(BIFIEobj, vars, group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.univar}(object,digits=3,...)
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
 \item{stat}{Data frame with univariate statistics}
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
See \code{\link{BIFIE.univar.test}} for a test of equal means and
effect sizes \eqn{\eta} and \eqn{d}.

\code{\link[survey:svymean]{svymean}} (\pkg{survey}),
\code{\link[intsvy:timss.mean]{timss.mean}} (\pkg{intsvy}),
\code{\link[intsvy:timss.mean.pv]{timss.mean.pv}} (\pkg{intsvy}),
\code{\link[stats:weighted.mean]{weighted.mean}} (\pkg{stats}),
\code{\link[Hmisc:wtd.stats]{wtd.mean}} (\pkg{Hmisc}),
\code{\link[miceadds:ma.wtd.statNA]{ma.wtd.meanNA}} (\pkg{miceadds})

\code{\link[survey:svyvar]{svyvar}} (\pkg{survey}),
\code{\link[Hmisc:wtd.stats]{wtd.var}} (\pkg{Hmisc}),
\code{\link[miceadds:ma.wtd.statNA]{ma.wtd.sdNA}} (\pkg{miceadds}),
\code{\link[miceadds:ma.wtd.statNA]{ma.wtd.covNA}} (\pkg{miceadds})
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

# compute descriptives for plausible values
res1 <- BIFIE.univar( bdat , vars = c("ASMMAT","ASSSCI","books") )
summary(res1)

# split descriptives by number of books
res2 <- BIFIE.univar( bdat , vars = c("ASMMAT","ASSSCI") , group="books" ,
            group_values=1:5)
summary(res2)

#############################################################################
# EXAMPLE 2: TIMSS dataset with missings
#############################################################################

data(data.timss2)
data(data.timssrep)

# use first dataset with missing data from data.timss2
bdat1 <- BIFIE.data( data.list=data.timss2[[1]] , wgt=  data.timss2[[1]]$TOTWGT ,
               wgtrep=data.timssrep[, -1 ])

# some descriptive statistics without statistical inference
res1a <- BIFIE.univar( bdat1 , vars = c("ASMMAT","ASSSCI","books") , se=FALSE)
# descriptive statistics with statistical inference
res1b <- BIFIE.univar( bdat1 , vars = c("ASMMAT","ASSSCI","books") )
summary(res1a)
summary(res1b)

# split descriptives by number of books
res2 <- BIFIE.univar( bdat1 , vars = c("ASMMAT","ASSSCI") , group="books" )
# Note that if group_values is not specified as an argument it will be
# automatically determined by the observed frequencies in the dataset 
summary(res2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}
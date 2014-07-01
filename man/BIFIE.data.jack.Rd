\name{BIFIE.data.jack}
\alias{BIFIE.data.jack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create \code{BIFIE.data} Object with Jackknife Zones
}

\description{
Creates a \code{BIFIE.data} objects for designs with jackknife zones,
especially for TIMSS/PIRLS studies.
}

\usage{
BIFIE.data.jack(data, wgt = NULL, pv_vars = NULL, jkzone = NULL, jkrep = NULL, 
     jktype = "JK_TIMSS", jkfac = NULL, fayfac = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame: Can be a single or a list of multiply imputed datasets
}
  \item{wgt}{
Case weight. If \code{jktype='JK_TIMSS'}, then \code{wgt='TOTWGT'}.
}
  \item{pv_vars}{
An optional vector of plausivle values which define multiply
imputed datasets.
}
  \item{jkzone}{
Jackknife zones. If \code{jktype='JK_TIMSS'}, then \code{jkzone='JKZONE'}.
}
  \item{jkrep}{
Jackknife replicate factors. If \code{jktype='JK_TIMSS'}, then \code{jkrep='JKREP'}.
}
  \item{jktype}{
Type of jackknife procdure for creating the \code{BIFIE.data} object.
\code{jktype='JK_TIMSS'} refers to TIMSS/PIRLS datasets.
}
  \item{jkfac}{
Factor for multiplying jackknife replicate weights. 
If \code{jktype='JK_TIMSS'}, then \code{jkfac=2}.
}

  \item{fayfac}{
Fay factor for statistical inference. The default is set to 1.
}
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
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BIFIE.data}}
}



\examples{
#############################################################################
# EXAMPLE 1: Convert TIMSS dataset to BIFIE.data object
#############################################################################

data(data.timss3)

# define plausible values
pv_vars <- c("ASMMAT" , "ASSSCI" )
# create BIFIE.data objects -> 5 imputed datasets
bdat1 <- BIFIE.data.jack( data=data.timss3  ,  pv_vars = pv_vars ,
             jktype="JK_TIMSS"  )
summary(bdat1)

# create BIFIE.data objects -> all PVs are included in one dataset
bdat2 <- BIFIE.data.jack( data=data.timss3  ,  jktype="JK_TIMSS"  )
summary(bdat2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIE data}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

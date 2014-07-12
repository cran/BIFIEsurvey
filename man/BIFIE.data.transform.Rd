\name{BIFIE.data.transform}
\alias{BIFIE.data.transform}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Data Transformation for \code{BIFIEdata} Objects
}

\description{
Computed a data transformation for \code{BIFIEdata} objects.
}

\usage{
BIFIE.data.transform( bifieobj , transform.formula )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bifieobj}{
Object of class \code{BIFIEdata}
}
  \item{transform.formula}{
\R formula object for data transformation
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
An object of class \code{BIFIEdata}. Additional values are
 \item{varnames.added}{Added variables in data transformation}
 \item{varsindex.added}{Indices of added variables}
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

%\seealso{
%\code{\link[survey:svyby]{svyby}} (\pkg{survey})
%}


\examples{
#############################################################################
# EXAMPLE 1: Data transformations for TIMSS data
#############################################################################

data(data.timss2)
data(data.timssrep)
# create BIFIEdata object
bifieobj <- BIFIE.data( data.timss2 , wgt=data.timss2[[1]]$TOTWGT , 
            wgtrep=data.timssrep[,-1] )

#*** Transformation 1: Squared book variable and interaction of migrant and female
transform.formula <- ~ 0 + I( books^2 ) + I( migrant * female )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula )
##   Included 2 variables: I(books^2) I(migrant * female)
# rename added variables
bifieobj$varnames[ bifieobj$varsindex.added ] <- c("books_sq" , "migrant_female")

\dontrun{
#*** Transformation 2: Create data variable for books
transform.formula <- ~ 0 + as.factor(books)
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula )
##   Included 5 variables: as.factor(books)1 as.factor(books)2 as.factor(books)3 
##        as.factor(books)4 as.factor(books)5
bifieobj$varnames[ bifieobj$varsindex.added ] <- paste0("books_D" , 1:5)

#*** Transformation 3: Discretized mathematics score
hi3a <- BIFIE.hist( bifieobj , vars = "ASMMAT" )
plot(hi3a)
transform.formula <- ~ 0 + I( as.numeric(cut( ASMMAT , breaks=seq(200,800, 100) )) )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula )
bifieobj$varnames[ bifieobj$varsindex.added ] <- "ASMMAT_discret"
hi3b <- BIFIE.hist( bifieobj , vars = "ASMMAT_discret" , breaks = 1:7 )
plot(hi3b)
# check frequencies
fr3b <- BIFIE.freq( bifieobj , vars = "ASMMAT_discret" , se=FALSE )
summary(fr3b)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data transformation}


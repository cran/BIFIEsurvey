\name{BIFIE.linreg}
\alias{BIFIE.linreg}
\alias{summary.BIFIE.linreg}
\alias{coef.BIFIE.linreg}
\alias{vcov.BIFIE.linreg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Linear Regression
}

\description{
Computes linear regression.
}

\usage{
BIFIE.linreg(BIFIEobj, dep=NULL , pre=NULL, formula=NULL , 
    group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.linreg}(object,digits=4,...)

\method{coef}{BIFIE.linreg}(object,...)

\method{vcov}{BIFIE.linreg}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
\item{dep}{String for the dependent variable in the regression model}
 \item{pre}{
Vector of predictor variables. If the intercept should be included,
then use the variable \code{one} for specifying it (see Examples).
}
\item{formula}{An \R formula object which can be applied instead of
	providing \code{dep} and \code{pre}. Note that there is 
   additional computation time needed for model matrix creation.}
  \item{group}{
Optional grouping variable
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{object}{Object of class \code{BIFIE.linreg}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with unstandardized and standardized regression
 coefficients, residual standard deviation and \eqn{R^2}}
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
\code{\link[survey:svyglm]{svyglm}} (\pkg{survey}),
\code{\link[intsvy:timss.reg]{timss.reg}} (\pkg{intsvy}),	
\code{\link[intsvy:timss.reg.pv]{timss.reg.pv}} (\pkg{intsvy}),
\code{\link[stats:lm]{lm}} (\pkg{stats})
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

#**** Model 1: Linear regression for mathematics score
mod1 <- BIFIE.linreg( bdat , dep= "ASMMAT" , pre=c("one","books","migrant") ,
              group= "female" )
summary(mod1)

\dontrun{
# same model but specified with R formulas 
mod1a <- BIFIE.linreg( bdat , formula = ASMMAT ~ books + migrant ,
               group= "female" , group_values = 0:1 )
summary(mod1a)

# compare result with lm function and first imputed dataset
dat1 <- data.timss1[[1]]
mod1b <- lm( ASMMAT ~ 0+as.factor(female)+as.factor(female):books+as.factor(female):migrant , 
                data= dat1 ,  weights=dat1$TOTWGT )
summary(mod1b)

#**** Model 2: Like Model 1, but books is now treated as a factor 
mod2 <- BIFIE.linreg( bdat , formula = ASMMAT ~ as.factor(books) + migrant  )
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linear regression}
\keyword{summary}
\keyword{coef}
\keyword{vcov}


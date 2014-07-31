\name{BIFIE.data.transform}
\alias{BIFIE.data.transform}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Data Transformation for \code{BIFIEdata} Objects
}

\description{
Computes a data transformation for \code{BIFIEdata} objects.
}

\usage{
BIFIE.data.transform( bifieobj , transform.formula , varnames.new=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bifieobj}{
Object of class \code{BIFIEdata}
}
  \item{transform.formula}{
\R formula object for data transformation.
}
\item{varnames.new}{Optional vector of names for new defined variables.
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
library(miceadds)

#############################################################################
# EXAMPLE 1: Data transformations for TIMSS data
#############################################################################

data(data.timss2)
data(data.timssrep)
# create BIFIEdata object
bifieobj1 <- BIFIE.data( data.timss2 , wgt=data.timss2[[1]]$TOTWGT , 
            wgtrep=data.timssrep[,-1] )
# create BIFIEdata object in compact way (cdata=TRUE)            
bifieobj2 <- BIFIE.data( data.timss2 , wgt=data.timss2[[1]]$TOTWGT , 
            wgtrep=data.timssrep[,-1] , cdata=TRUE)
            
#*** Transformation 1: Squared and cubic book variable
transform.formula <- ~ 0 + I( books^2 ) + I( books^3 )
# as.character(transform.formula)
bifieobj <- BIFIE.data.transform( bifieobj1 , transform.formula =transform.formula )
bifieobj$variables
# rename added variables
bifieobj$varnames[ bifieobj$varsindex.added ] <- c("books_sq" , "books_cub")

# check descriptive statistics
res1 <- BIFIE.univar( bifieobj , vars = c("books_sq" , "books_cub" ) )
summary(res1)

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
transform.formula <- ~ 0 + I( as.numeric(cut( ASMMAT , breaks=seq(200,800,100) )) )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula ,
                  varnames.new = "ASMMAT_discret" )
hi3b <- BIFIE.hist( bifieobj , vars = "ASMMAT_discret" , breaks = 1:7 )
plot(hi3b)
# check frequencies
fr3b <- BIFIE.freq( bifieobj , vars = "ASMMAT_discret" , se=FALSE )
summary(fr3b)

#*** Transformation 4: include standardization variables for book variable

# start with testing the transformation function on a single dataset
dat1 <- bifieobj$dat1
weighted.mean( dat1[,"books"] , dat1[,"TOTWGT"] , na.rm=TRUE)
sqrt( Hmisc::wtd.var( dat1[,"books"] , dat1[,"TOTWGT"] , na.rm=TRUE) )
# z standardization
transform.formula <- ~ 0 + I( ( books - weighted.mean( books , TOTWGT , na.rm=TRUE) ) /
                                sqrt( Hmisc::wtd.var( books , TOTWGT , na.rm=TRUE) )  )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula ,
                  varnames.new = "z_books" )
# standardize variable books with M=500 and SD=100
transform.formula <- ~ 0 + I( 500 + 100*( books - weighted.mean( books , TOTWGT , na.rm=TRUE) ) /
                                sqrt( Hmisc::wtd.var( books , TOTWGT , na.rm=TRUE) )  )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula ,
                  varnames.new = "z500_books" )
                  
# standardize variable books with respect to M and SD of ALL imputed datasets
res <- BIFIE.univar( bifieobj , vars = "books" )
summary(res)
##       var  Nweight Ncases     M M_SE M_fmi M_VarMI M_VarRep    SD SD_SE SD_fmi
##   1 books 76588.72   4554 2.945 0.04     0       0    0.002 1.146 0.015      0
M <- round(res$output$mean1,5)
SD <- round(res$output$sd1,5)
transform.formula <- paste0( " ~ 0 + I( ( books - " ,  M  , " ) / " , SD  , ")"  )
##   > transform.formula
##   [1] " ~ 0 + I( ( books - 2.94496 ) / 1.14609)"
bifieobj <- BIFIE.data.transform( bifieobj, transform.formula =as.formula(transform.formula),
                 varnames.new = "zall_books" )                  
                  
# check statistics
res4 <- BIFIE.univar( bifieobj , vars = c("z_books" , "z500_books" , "zall_books") )
summary(res4)

#*** Transformation 5: include rank transformation for variable ASMMAT

# calculate percentage ranks using wtd.rank function from Hmisc package
dat1 <- bifieobj$dat1
100 * Hmisc::wtd.rank( dat1[,"ASMMAT"] , w=dat1[,"TOTWGT"] ) / sum( dat1[,"TOTWGT"] )
# define an auxiliary function for calculating percentage ranks
wtd.percrank <- function( x , w ){
        100 * Hmisc::wtd.rank( x , w , na.rm=TRUE ) / sum( w , na.rm=TRUE )
                    }
wtd.percrank( dat1[,"ASMMAT"] , dat1[,"TOTWGT"] )
# define transformation formula
transform.formula <- ~ 0 + I( wtd.percrank( ASMMAT , TOTWGT ) )
# add ranks to BIFIEdata object
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula ,
                  varnames.new = "ASMMAT_rk" )
# check statistic
res5 <- BIFIE.univar( bifieobj , vars = c("ASMMAT_rk" ) )
summary(res5)

#*** Transformation 6: recode variable books

# recode variable books according to "1,2=0 , 3,4=1 , 5=2"
dat1 <- bifieobj$dat1
# use Recode function from car package
car::Recode( dat1[,"books"], "1:2='0'; c(3,4)='1';5='2'")
# define transformation formula
transform.formula <- ~ 0 + I( car::Recode( books, "1:2='0'; c(3,4)='1';5='2'") )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula =transform.formula ,
                  varnames.new = "book_rec" )
res6 <- BIFIE.freq( bifieobj , vars = c("book_rec" ) )
summary(res6)

#*** Transformation 7: include some variables aggregated to the school level

dat1 <- as.data.frame(bifieobj$dat1)
# at first, create school ID in the dataset by transforming the student ID
dat1$idschool <- as.numeric(substring( dat1$IDSTUD , 1 , 5 ))
transform.formula <- ~ 0 +I( as.numeric( substring( IDSTUD  , 1 , 5 ) ) )
bifieobj <- BIFIE.data.transform( bifieobj , transform.formula = transform.formula  ,
               varnames.new = "idschool" )
                          
# define group mean for variable ASMMAT
dat1 <- as.data.frame(bifieobj$dat1)
gm <- miceadds::fast.groupmean( data=dat1$ASMMAT , group=dat1$idschool )
gm[ match( dat1$idschool , gm[,1] ) , 2 ]   # school means replicated for each student

# define utility function for group mean
add.groupmean <- function( x , group ){
    gm <- miceadds::fast.groupmean( data=x , group=group )
    gm[ match( group , gm[,1] ) , 2 ]  
                }
# add school mean ASMMAT
bifieobj <- BIFIE.data.transform( bifieobj , 
                transform.formula = ~ 0 + I( add.groupmean( ASMMAT , idschool ) ) , 
                varnames.new = "M_ASMMAT" )
# add school mean books
bifieobj <- BIFIE.data.transform( bifieobj , 
                transform.formula = ~ 0 + I( add.groupmean( books , idschool ) ) , 
                varnames.new = "M_books" )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data transformation}


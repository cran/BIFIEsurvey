\name{BIFIEsurvey-package}
\alias{BIFIEsurvey-package}
\alias{BIFIEsurvey}
\docType{package}

\title{
Some Tools for Survey Statistics in Educational Assessment, Developed by BIFIE
}
\description{
    Some Tools for Survey Statistics in Educational Assessment, 
    Developed by BIFIE (\url{http://www.bifie.at}).
    Bundesinstitut fuer Bildungsforschung,
    Entwicklung und Innovation im oesterreichischen Schulwesen;
    Federal Institute for Educational Research, Innovation and 
    Development of the Austrian School System;
    Salzburg, Austria.
}
\details{
\tabular{ll}{
Package: \tab \pkg{BIFIEsurvey} \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014 \cr
License: \tab GPL (>=2)\cr
}

This package include basic descriptive functions for large scale assessment data 
to complement the more comprehensive \pkg{survey} package. The few
functions in this package were written in \pkg{Rcpp}.

The features of \pkg{BIFIEsurvey} include for designs with replicate weights 
(which includes Jackknife and Bootstrap as general approaches):

\itemize{

\item Descriptive statistics: means and standard deviations (\code{\link{BIFIE.univar}}),
frequencies (\code{\link{BIFIE.freq}}), 
crosstabs (\code{\link{BIFIE.crosstab}})

\item Linear regression (\code{\link{BIFIE.linreg}})

\item Logistic regression (\code{\link{BIFIE.logistreg}})

\item Statistical inference for derived parameters (\code{\link{BIFIE.derivedParameters}})

\item Wald tests (\code{\link{BIFIE.waldtest}}) of model parameters based on 
replicated statistics


%% logistic regression
%% \item jackknife and bootstrap for non-survey data
	
}




}



\author{
BIFIE
\url{http://www.bifie.at}

Maintainer: Alexander Robitzsch <a.robitzsch@bifie.at>
}

%----
%For the following functions, an own function maintainer is defined
%which should be contacted instead of the package maintainer: 
%
%Alexander Robitzsch <a.robitzsch@bifie.at> \cr
%\code{\link{BIFIE.logistreg}}
%----

%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{package}
\seealso{
See also the \pkg{survey}, \pkg{intsvy}, \pkg{svyPVpack},
\pkg{lavaan.survey}, \pkg{EVER} and 
the \pkg{eatRep} (\url{https://r-forge.r-project.org/R/?group_id=1326})
packages.
}

\examples{
##   |-----------------------------------------------------------------
##   | BIFIEsurvey 0.1-21 (2014-06-21)                                       
##   | Maintainer: Alexander Robitzsch <a.robitzsch at bifie.at >      
##   | http://www.bifie.at                                             
##   |-----------------------------------------------------------------
  
##   .........................*  ,::;              :,:;       * .;*;. .,:                
##   :::::::::::::::::::::::::.  ##+@              ##+#       .@####+ ;+# *              
##   :::::::::::::::::::::::::.  ###@              ####       @@;   :*;##**              
##   :::::::::::::::::::::::::.  ###@              ##+#     *##.      .,,                
##   :::::::::::::::::::::::::.  ###@              ::,:    * ##                          
##   :::::::::::::::::::::::::.  ###@                 *      @@                          
##   :::::::::::::::::::::::::.  ###@                      * #@                          
##   :::::::::::::::::::::::::.  ###@                        #@                   *      
##   :::::::::::::::::::::::::.  ##@# ,@###@       @@## * @@@+##@@@@   #@     *.@#####   
##   :::::::::::::::::::::::::.  ####*@#####@.**   ####   @@#@@@####   #@     ;##+**#+@* 
##   :::::::::::::::::::::::::.  ##@@##,,#+##@     ####      @@        #@    .#@*  * *## 
##   ::::,::::::::::::::::::::.  ##+@@    ####,    ####      #@        #@    ##,      #+*
##   ::::**.::::::::::::::::::.  ##+@, ** ,###@    ####      #@        #@    @@       ;@;
##   :::::*  .::::::::::::::::.  ##*       ####    ####      #@        #@   .##@@@@@@@#+#
##   :::::::*   *,,:::::::::::.     :      ###@    ####      #@        #@   ;##@@@@@@@@@@
##   :::::::::.      **....*      ,@#      ####    ####      #@        #@   *@**         
##   ::::::::::::.*              ####      @###*   ####      #@        #@  *+#*          
##   ::::::::::::::::,,,,,,,::.  ####      ####    ####      #@        #@  **#*          
##   :::::::::::::::::::::::::.  ####      @###    ####      #@        #@   ,##          
##   :::::::::::::::::::::::::.  ###@   *  @##+    ####      #@        #@   *#@          
##   :::::::::::::::::::::::::*  @@##     ,#@##    ####      #@        #@    @#          
##   :::::::::::::::::::::::::*  @+##,    @###*    ####      #@        #@   *+#@         
##   ::::::::::::::::::::::::.    ###@:  *###@     ####      #@        #@     #@+        
##   :::::::::::::::::::::::.   **;@#@#@####@.     ####      #@        #@     *@#@:*   * 
##   ::::::::::::::::::::::.       ,@######@.      ####*     ##        @+      *#####@## 
##   ::::::::::::::::::::.*        * .*##*. *        ***     *.        **         ;##+;. 

}

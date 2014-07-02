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
Package: \tab BIFIEsurvey\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014 \cr
License: \tab GPL (>=2)\cr
}

This package include basic descriptive functions for large scale assessment data 
to complement the more comprehensive \pkg{survey} package. The few
functions in this package were written in \pkg{Rcpp}.

The features of \pkg{BIFIEsurvey} include for designs with replicate weights:

\itemize{

\item Descriptive statistics: means, standard deviations, frequencies,
crosstabs and Wald tests of these parameters based on replicated statistics

\item Linear regression

%% logistic regression
%% \item jackknife and bootstrap for non-survey data
	
}




}



\author{
BIFIE
\url{http://www.bifie.at}

Maintainer: Alexander Robitzsch <a.robitzsch@bifie.at>
}
%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{package}
\seealso{
See also the \pkg{survey}, \pkg{intsvy}, \pkg{svyPVpack},
\pkg{lavaan.survey} and 
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
\section{\R Function Versions }{
\emph{BIFIE.by__0.02.R}, 
\emph{BIFIE.correl__0.04.R}, 
\emph{BIFIE.crosstab__0.05.R}, \cr
\emph{BIFIE.data.jack__1.04.R}, 
\emph{BIFIE.data__1.06.R}, 
\emph{BIFIE.freq__0.08.R}, 
\emph{BIFIE.linreg__0.12.R}, \cr
\emph{BIFIE.progressbar__0.02.R}, 
\emph{BIFIE.summary__0.03.R}, 
\emph{BIFIE.univar.test__0.10.R}, \cr
\emph{BIFIE.univar__1.32.R}, 
\emph{BIFIE.waldtest__0.22.R}, 
\emph{print.object.summary__0.02.R}, \cr
\emph{zzz__0.02.R}, 
}
\section{\emph{Rcpp} Function Versions }{
\emph{jack_dataprep__0.02.cpp}, 
\emph{univar__3.23.cpp}, 
\emph{univar_helpers__2.70.h}, 
}
\section{\emph{Rd} Documentation Versions }{
\emph{BIFIE.by__0.06.Rd}, 
\emph{BIFIE.correl__0.09.Rd}, 
\emph{BIFIE.crosstab__0.04.Rd}, \cr
\emph{BIFIE.data.jack__0.09.Rd}, 
\emph{BIFIE.data__0.09.Rd}, 
\emph{BIFIE.freq__0.06.Rd}, \cr
\emph{BIFIE.linreg__0.12.Rd}, 
\emph{BIFIE.univar.test__0.05.Rd}, 
\emph{BIFIE.univar__0.18.Rd}, \cr
\emph{BIFIE.waldtest__0.16.Rd}, 
\emph{BIFIEsurvey-package__0.25.Rd}, 
\emph{data.timss1__0.08.Rd}, 
}

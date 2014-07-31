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
Version: \tab 0.3\cr
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
crosstabs (\code{\link{BIFIE.crosstab}}) and Wald tests (\code{\link{BIFIE.waldtest}})
of these parameters based on replicated statistics

\item Linear regression (\code{\link{BIFIE.linreg}})

\item Logistic regression (\code{\link{BIFIE.logistreg}})

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
\section{\R Function Versions }{
\emph{BIFIE.bifiedata2bifiecdata__2.04.R}, 
\emph{BIFIE.by__1.05.R}, 
\emph{BIFIE.correl__0.10.R}, \cr
\emph{BIFIE.crosstab__0.11.R}, 
\emph{BIFIE.data.boot__0.05.R}, 
\emph{BIFIE.data.jack__1.21.R}, \cr
\emph{BIFIE.data.select__1.05.R}, 
\emph{BIFIE.data.transform__1.08.R}, 
\emph{BIFIE.data__1.26.R}, \cr
\emph{BIFIE.ecdf__0.11.R}, 
\emph{BIFIE.freq__0.13.R}, 
\emph{BIFIE.hist__0.08.R}, 
\emph{BIFIE.linreg__0.17.R}, \cr
\emph{BIFIE.logistreg__0.15.R}, 
\emph{BIFIE.mva__0.05.R}, 
\emph{BIFIE.progressbar__0.02.R}, \cr
\emph{BIFIE.summary__0.03.R}, 
\emph{BIFIE.univar.test__0.10.R}, 
\emph{BIFIE.univar__1.37.R}, \cr
\emph{BIFIE.waldtest__0.23.R}, 
\emph{cdata.wgtrep__0.01.R}, 
\emph{load.BIFIE.data__1.14.R}, \cr
\emph{load.BIFIEdata.files__3.01.R}, 
\emph{print.object.summary__0.02.R}, 
\emph{save.BIFIE.data__1.01.R}, \cr
\emph{write.BIFIE.data__1.09.R}, 
\emph{zzz__0.02.R}, 
}
\section{\emph{Rcpp} Function Versions }{
\emph{jack_dataprep__0.23.cpp}, 
\emph{univar__3.57.cpp}, 
\emph{univar_helpers__3.35.h}, 
}
\section{\emph{Rd} Documentation Versions }{
\emph{BIFIE.BIFIEdata2BIFIEcdata__0.03.Rd}, 
\emph{BIFIE.by__0.07.Rd}, 
\emph{BIFIE.correl__0.10.Rd}, \cr
\emph{BIFIE.crosstab__0.04.Rd}, 
\emph{BIFIE.data.boot__0.06.Rd}, 
\emph{BIFIE.data.jack__0.19.Rd}, \cr
\emph{BIFIE.data.select__0.05.Rd}, 
\emph{BIFIE.data.transform__0.16.Rd}, 
\emph{BIFIE.data__0.22.Rd}, \cr
\emph{BIFIE.ecdf__0.11.Rd}, 
\emph{BIFIE.freq__0.07.Rd}, 
\emph{BIFIE.hist__0.13.Rd}, 
\emph{BIFIE.linreg__0.13.Rd}, \cr
\emph{BIFIE.logistreg__0.09.Rd}, 
\emph{BIFIE.mva__0.02.Rd}, 
\emph{BIFIE.univar.test__0.06.Rd}, \cr
\emph{BIFIE.univar__0.18.Rd}, 
\emph{BIFIE.waldtest__0.18.Rd}, 
\emph{BIFIEsurvey-package__0.38.Rd}, \cr
\emph{data.pisaNLD__0.07.Rd}, 
\emph{data.test1__0.02.Rd}, 
\emph{data.timss1__0.11.Rd}, \cr
\emph{save.BIFIEdata__0.30.Rd}, 
}

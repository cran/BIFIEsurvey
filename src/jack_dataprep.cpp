
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;



// declarations
extern "C" {
SEXP bifie_jack_timss( SEXP wgt_, SEXP jkzone_, SEXP jkrep_, SEXP RR_, 
	  SEXP jkfac_, SEXP prbar_) ;
}

// definition

SEXP bifie_jack_timss( SEXP wgt_, SEXP jkzone_, SEXP jkrep_, SEXP RR_, 
	 SEXP jkfac_, SEXP prbar_ ){
BEGIN_RCPP
   
       
       
     Rcpp::NumericVector wgt(wgt_);  
     Rcpp::NumericVector jkzone(jkzone_);  
     Rcpp::NumericVector jkrep(jkrep_);  
     int RR = as<int>(RR_);  
     double jkfac = as<double>(jkfac_) ;  
     Rcpp::NumericVector prbar(prbar_);  
       
     int N = wgt.size() ;  
       
     Rcpp::NumericMatrix wgtrep(N,RR);  
       
     for (int rr=0;rr<RR;rr++){  
     for (int nn=0;nn<N;nn++){  
     if ( jkzone[nn] == rr ){  
     	wgtrep(nn,rr) = jkfac * wgt[nn] * jkrep[nn] ;  
     			} else {  
     	wgtrep(nn,rr) = wgt[nn] ;  
     			}  
     		}  
     if (prbar[rr] == 1){		  
     	Rcpp::Rcout << "-" <<  std::flush  ;  
     		}  
     		  
     	}  
       
     return( wrap( wgtrep) ) ;  
       
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;
END_RCPP
}




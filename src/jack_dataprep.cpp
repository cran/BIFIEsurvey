
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



//*************************************************************************
//*************************************************************************
//*************************************************************************
//  generate Jackknife replicate weights TIMSS studies
//*************************************************************************
//*************************************************************************
//*************************************************************************



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

//*************************************************************************
//*************************************************************************
//*************************************************************************
// generate bootstrap samples (IID assumption)
//*************************************************************************
//*************************************************************************
//*************************************************************************


// declarations
extern "C" {
SEXP bifie_boot( SEXP cumwgt_, SEXP rand_wgt_) ;
}

// definition

SEXP bifie_boot( SEXP cumwgt_, SEXP rand_wgt_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericVector cumwgt(cumwgt_);          
     Rcpp::NumericMatrix rand_wgt(rand_wgt_) ;  
       
     int N = cumwgt.size() ;  
     int WW = rand_wgt.ncol() ;  
       
     arma::colvec wgt_temp(N);  
     Rcpp::NumericMatrix wgtM(N,WW);  
       
     int zz=0;  
     int nn1=0;  
       
     for (int ww=0;ww<WW;ww++){  
     for (int nn=0;nn<N;nn++){  
     	wgt_temp(nn,0) = rand_wgt(nn,ww) ;  
     			}  
     // sort entries  
     arma::colvec wgt_tempsort = arma::sort( wgt_temp ) ;  
       
     zz=0;  
     nn1=0;  
       
     while (nn1<N){  
     	if (wgt_tempsort(nn1,0) < cumwgt[zz] ){  
     		wgtM(zz,ww) ++ ;  
     		nn1 ++ ;	  
     		} else {  
     		zz ++ ;   
     		}  
     	}  
     }  
     	  
     //*************************************************      
     // OUTPUT                                     
     return Rcpp::List::create(   
         Rcpp::_["cumwgt"] = cumwgt  ,  
         Rcpp::_["rand_wgt"] = rand_wgt ,  
         Rcpp::_["wgtM"] = wgtM      
         ) ;    
     // maximal list length is 20!  
              
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
                  
END_RCPP
}

//*************************************************************************
//*************************************************************************
//*************************************************************************
// converting BIFIEdata object into BIFIEcdata object
//*************************************************************************
//*************************************************************************
//*************************************************************************


// declarations
extern "C" {
SEXP bifie_bifiedata2bifiecdata( SEXP datalistM_, SEXP Nimp_) ;
}

// definition

SEXP bifie_bifiedata2bifiecdata( SEXP datalistM_, SEXP Nimp_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalistM(datalistM_);          
     int Nimp=as<int>(Nimp_);       
       
     int N= datalistM.nrow() / Nimp ;  
     int VV = datalistM.ncol() ;  
       
     //************************************************  
     //*** define indicator matrix  
     Rcpp::NumericMatrix datalistM_ind(N, VV);  
     double t1=0;  
     // number of missing entries  
     int Nmiss=0;  
     for (int nn=0;nn<N;nn++){ // beg nn  
     for (int vv=0;vv<VV;vv++){  // beg vv  
     	datalistM_ind(nn,vv) = 1 ;  
     	t1 = datalistM(nn,vv) ;  
     if ( ! R_IsNA( datalistM(nn , vv ) ) ){		  
     	for (int ii=1;ii<Nimp;ii++){ // beg ii  
     		if ( datalistM( nn + ii*N , vv ) != t1 ){  
     			datalistM_ind(nn,vv)=0;  
     			Nmiss ++ ;  
     			break ;  
     				}  
     			}  // end ii  
     		}  
     	} // end vv	  
     	} // end nn  
       
       
     //************************************************  
     // imputed data matrix  
     // int ZZ=Nmiss*Nimp;  
     // Rcpp::NumericMatrix datalistM_imputed(ZZ,4);  
     // first column: imputed dataset  
     // second column: subject case  
     // third column: imputed variable  
     // fourth column: value  
       
     // restructuring ARb 2014-07-30  
     int ZZ=Nmiss;  
     Rcpp::NumericMatrix datalistM_imputed(ZZ,Nimp);  
     Rcpp::IntegerMatrix datalistM_impindex(ZZ,2);  
     // first column: subject case  
     // second column: imputed variable  
       
       
     int zz=0;  
     for (int vv=0;vv<VV;vv++){  
     for (int nn=0;nn<N;nn++){  
     if ( ! R_IsNA( datalistM(nn , vv ) ) ){	  
     if ( datalistM_ind(nn,vv) == 0 ){ // beg if missing entry  
     	datalistM_impindex(zz,0)=nn;  
     	datalistM_impindex(zz,1)=vv;				  
     	for (int ii=0;ii<Nimp;ii++){  // beg ii  
     	//	datalistM_imputed(zz,0) = ii;  
     	//	datalistM_imputed(zz,1) = nn;	  
     	//	datalistM_imputed(zz,2) = vv;	  
     	//	datalistM_imputed(zz,3) = datalistM( nn+ii*N , vv) ;  
     	datalistM_imputed(zz,ii) = datalistM(nn+ii*N , vv ) ;		  
     				} // end ii  
     	zz ++ ;			  
     		}	// end if missing entry  
     }  
     } // end nn  
     } // end vv  
       
     	  
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         Rcpp::_["datalistM_ind"] = datalistM_ind ,  
         Rcpp::_["datalistM_imputed"] = datalistM_imputed ,  
         Rcpp::_["datalistM_impindex"] = datalistM_impindex ,      
         Rcpp::_["Nimp"] = Nimp ,  
         Rcpp::_["Nmiss"] = Nmiss  
         ) ;    
       
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}


//*************************************************************************
//*************************************************************************
//*************************************************************************
// convert BIFIEcdata to BIFIEdata
//*************************************************************************
//*************************************************************************
//*************************************************************************


// declarations
extern "C" {
SEXP bifie_bifiecdata2bifiedata( SEXP datalistM_ind_, SEXP datalistM_imputed_, 
	SEXP Nimp_, SEXP dat1_, SEXP datalistM_impindex_) ;
}

// definition

SEXP bifie_bifiecdata2bifiedata( SEXP datalistM_ind_, SEXP datalistM_imputed_, 
	SEXP Nimp_, SEXP dat1_, SEXP datalistM_impindex_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalistM_ind(datalistM_ind_);  
     Rcpp::NumericMatrix datalistM_imputed(datalistM_imputed_);  
     int Nimp=as<int>(Nimp_);       
     Rcpp::NumericMatrix dat1(dat1_);  
     Rcpp::NumericMatrix datalistM_impindex(datalistM_impindex_);  
       
       
     int N=dat1.nrow();  
     int VV=dat1.ncol();  
     int ZZ=Nimp*N;  
       
     Rcpp::NumericMatrix datalistM(ZZ,VV);  
       
     //**** non-imputed data  
     for (int ii=0;ii<Nimp;ii++){   // beg ss  
     for (int nn=0;nn<N;nn++){ // beg nn  
        for (int vv=0;vv<VV;vv++){	// beg vv  
     	if ( datalistM_ind(nn,vv) == 1 ){ // beg non-imputed data  
     		datalistM(nn+ii*N,vv) = dat1(nn,vv) ;  
     				}  // end non-imputed data			  
     			} // end vv  
     		} // end nn  
           } // end ss  
       
             
     //**** imputed data  
     //--- Rcpp::NumericMatrix datalistM_imputed(ZZ,4);  
     // 	first column: subject case  
     // 	second column: imputed variable  
       
       
     int HH=datalistM_imputed.nrow();  
     // int ii_ ;  
     int nn_ ;  
     int vv_ ;  
       
     for ( int hh=0;hh<HH;hh++){  
        nn_ = datalistM_impindex(hh,0) ;	  
        vv_ = datalistM_impindex(hh,1) ;	  
        for (int ii=0;ii<Nimp;ii++){		  
     //	ii_ = datalistM_imputed(hh,0) ;  
             datalistM( nn_ + ii * N , vv_ ) = datalistM_imputed(hh, ii ) ;  
     			}	  
     }                                  
             
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         Rcpp::_["datalistM"] = datalistM  ,       
         Rcpp::_["Nimp"] = Nimp     
         ) ;    
       
     // maximal list length is 20!         
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;                                                        
     
END_RCPP
}



//*************************************************************************
//*************************************************************************
//*************************************************************************
//******************************************************************
// data preparation loading files with indicator matrix
//*************************************************************************
//*************************************************************************
//*************************************************************************


// declarations
extern "C" {
SEXP bifie_bifiedata_stepwise( SEXP dat_, SEXP dat_ind_, SEXP Nmiss_) ;
}

// definition

SEXP bifie_bifiedata_stepwise( SEXP dat_, SEXP dat_ind_, SEXP Nmiss_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix dat1(dat_);  
     Rcpp::NumericMatrix dat_ind(dat_ind_);  
     int Nmiss = as<int>(Nmiss_) ;  
       
     int N=dat1.nrow();  
     int VV=dat1.ncol();  
       
     //************************************************  
     // imputed data matrix  
     int ZZ=Nmiss;  
     Rcpp::NumericMatrix datalistM_imputed(ZZ,4);  
     // first column: imputed dataset  
     // second column: subject case  
     // third column: imputed variable  
     // fourth column: value  
     int zz=0;  
     int ii=0;  
     for (int vv=0;vv<VV;vv++){  
     for (int nn=0;nn<N;nn++){  
     if ( dat_ind(nn,vv) == 0 ){ // beg if missing entry  
     	datalistM_imputed(zz,0) = ii;  
     	datalistM_imputed(zz,1) = nn;	  
     	datalistM_imputed(zz,2) = vv;	  
     	datalistM_imputed(zz,3) = dat1( nn , vv) ;  
     	zz ++ ;				  
     	}	// end if missing entry  
     } // end nn  
     } // end vv  
       
       
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         Rcpp::_["datalistM_imputed"] = datalistM_imputed // ,      
         ) ;    
       
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}










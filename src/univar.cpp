


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

#include "univar_helpers.h"
// user includes



// declarations
extern "C" {
SEXP univar_multiple_V2group( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, 
	  SEXP vars_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_) ;
}

// definition

SEXP univar_multiple_V2group( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, 
	  SEXP vars_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
       
     int Nimp = NI[0] ;  
     int RR = wgtrep.ncol() ;   
     int N = wgt1.nrow() ;  
     int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
     int VV1=VV*GG ;  
       
     Rcpp::NumericMatrix mean1M(VV1,Nimp);  
     Rcpp::NumericMatrix sd1M(VV1,Nimp);  
     Rcpp::NumericMatrix mean1_varM(VV1,Nimp);  
     Rcpp::NumericMatrix sd1_varM(VV1,Nimp);  
     Rcpp::NumericMatrix mean1repM(VV1,RR*Nimp);  
     Rcpp::NumericMatrix sd1repM(VV1,RR*Nimp);  
     Rcpp::NumericMatrix dat1(N,NV) ;  
     Rcpp::NumericVector sumweights(1);  
     Rcpp::NumericVector mean1(VV1);  
     Rcpp::NumericVector sd1(VV1);  
     // int NP=mean1M.nrow();  
     Rcpp::NumericVector sumwgt1(VV);  
     Rcpp::NumericVector ncases1(VV*GG);  
     Rcpp::NumericVector mean1_var(VV1);   
     Rcpp::NumericVector sd1_var(VV1);  
     Rcpp::NumericMatrix sumweightM(GG,Nimp) ;  
     Rcpp::NumericMatrix ncasesM(GG*VV,Nimp) ;  
     Rcpp::NumericMatrix mean1rep(VV1,RR);  
     Rcpp::NumericMatrix sd1rep(VV1,RR);  
     Rcpp::NumericMatrix sumweightrepM(GG,RR*Nimp) ;
     
     Rcpp::Rcout << "|"  ;  
       
     //***********************  
     // loop multiply imputed datasets  
     for ( int ii=0; ii < Nimp ; ii++){  
       
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
     // statistics single data  
     Rcpp::List res1 = univar_helper_multiple_V2group( dat1,  wgt1,  vars_index ,  
     	group_index1 , group_values ) ;  
     mean1 = matr2vec(res1["mean1"]) ;  
     sd1 = matr2vec(res1["sd1"]) ;  
     sumwgt1 = matr2vec(res1["sumwgt1"]) ;  
     ncases1 = matr2vec(res1["ncases1"]) ;  
       
     // compute statistics for replicate weights ;  
     Rcpp::List res2 = univar_helper_multiple_V2group( dat1,  wgtrep,  vars_index ,  
     	group_index1 , group_values ) ;  
     Rcpp::NumericMatrix mean1rep = res2["mean1"] ;  
     Rcpp::NumericMatrix sd1rep = res2["sd1"] ;  
//     Rcpp::NumericMatrix sd1rep = res2["sd1"] ;
     Rcpp::NumericMatrix sumwgt1a = res2["sumwgt1"] ;     
     
         // compute standard errors  
         mean1_var = varjack_helper( mean1 , mean1rep , fayfac ) ;  
         sd1_var = varjack_helper( sd1 , sd1rep , fayfac ) ;		  
       
         // collect all results for one imputed dataset  
         for (int vv=0;vv<VV1;vv++){  
            mean1M( vv , ii ) = mean1[vv] ;  
            sd1M( vv , ii ) = sd1[vv] ;  
            mean1_varM( vv , ii ) = mean1_var[vv] ;  
            sd1_varM( vv , ii ) = sd1_var[vv] ;   
            for (int rr=0;rr<RR;rr++){  
     	    mean1repM( vv , rr+ii*RR ) = mean1rep(vv,rr) ;    
     	    sd1repM( vv , rr+ii*RR ) = sd1rep(vv,rr) ;  
     	   			   }  
     			}  
       
            for (int gg=0;gg<GG;gg++){  
            	      sumweightM(gg,ii) = sumwgt1[gg] ;
            	      for (int rr=0;rr<RR;rr++){
            	      	sumweightrepM( gg , rr+ii*RR) = sumwgt1a(gg,rr) ;      
            	      }
            	      for (int vv=0;vv<VV;vv++){
            	         ncasesM(gg+vv*GG,ii) = ncases1[gg+vv*GG] ;
            	         		}
            	       			}			  
       
        Rcpp::Rcout << "-" <<  std::flush ;   
        	// << std::endl ;       	       			  
            	       			  
     	} // end loop ii | imputed datasets	  
       
        Rcpp::Rcout << "|" << std::endl ;  	  
     	  
     //----  
     // inference multiply imputed datasets  
       
     Rcpp::List res3 = rubin_rules_univ( mean1M , mean1_varM ) ;  
     mean1=res3["pars"] ;  
     Rcpp::NumericVector mean1_se=res3["pars_se"] ;  
     Rcpp::NumericVector mean1_varWithin=res3["pars_varWithin"] ;  
     Rcpp::NumericVector mean1_varBetween=res3["pars_varBetween"] ;  
     Rcpp::NumericVector mean1_fmi=res3["pars_fmi"] ;  
       
     res3 = rubin_rules_univ( sd1M , sd1_varM ) ;  
     sd1=res3["pars"] ;  
     Rcpp::NumericVector sd1_se=res3["pars_se"] ;  
     Rcpp::NumericVector sd1_varWithin=res3["pars_varWithin"] ;  
     Rcpp::NumericVector sd1_varBetween=res3["pars_varBetween"] ;  
     Rcpp::NumericVector sd1_fmi=res3["pars_fmi"] ;	  
     	  
//     res3 = rubin_rules_univ( sumweightM , sumweightM ) ;  
//     Rcpp::NumericVector sumweight =res3["pars"] ;  
       
     res3 = rubin_rules_univ( ncasesM , ncasesM ) ;  
     Rcpp::NumericVector ncases =res3["pars"] ;  
       
     // Rcpp::Rcout << "m1repM " <<  mean1repM(3,10) <<  std::endl ;  
       
     //*************************************************      
     // OUTPUT              
       
                
     return Rcpp::List::create(   
         _["mean1"] = mean1 ,  
         _["mean1_se"] = mean1_se ,  
         _["mean1_varWithin"] = mean1_varWithin ,  
         _["mean1_varBetween"] = mean1_varBetween ,  
         _["mean1_fmi"] = mean1_fmi ,      
         _["mean1M"] = mean1M ,   
         _["mean1_varM"] = mean1_varM   ,  
         _["mean1repM"] = mean1repM ,  
         _["sd1"] = sd1 ,  
         _["sd1_se"] = sd1_se ,       
         _["sd1_varWithin"] = sd1_varWithin ,  
         _["sd1_varBetween"] = sd1_varBetween ,  
         _["sd1_fmi"] = sd1_fmi ,      
         _["sd1M"] = sd1M ,   _["sd1_varM"] = sd1_varM  ,   _["sd1repM"] = sd1repM ,  
         _["sumweightrepM"] = sumweightrepM ,  _["sumweightM"] = sumweightM  ,      
         _["ncases"] = ncases  ,  _["ncasesM"] = ncasesM    
         ) ;    
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;
END_RCPP
}






// declarations
extern "C" {
SEXP bifie_freq( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
   SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP vars_values_, SEXP vars_values_numb_) ;
}

// definition

SEXP bifie_freq( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP vars_values_, SEXP vars_values_numb_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
     Rcpp::NumericMatrix vars_values(vars_values_) ;  
     Rcpp::NumericVector vars_values_numb(vars_values_numb_) ;  
       
       
     int Nimp = NI[0] ;  
     int RR = wgtrep.ncol() ;   
     int N = wgt1.nrow() ;  
     int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     // number of values per variable  
     int VV1=0;  
     Rcpp::NumericVector vars_values_numb_cumsum(VV+1) ;  
     for (int vv=0;vv<VV; vv++){  
        VV1 += vars_values_numb[vv] ;  
        vars_values_numb_cumsum[vv+1] = vars_values_numb_cumsum[vv] + vars_values_numb[vv] ;  
        		}  
     int VV2=VV1*GG ;  
       
     // matrices for output  
     int WW = wgt1.ncol();   
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
     Rcpp::NumericMatrix perc1(VV2,WW) ;  
     Rcpp::NumericMatrix perc2(VV2,WW) ;  
     Rcpp::NumericVector sumwgt(VV*GG) ;  
     Rcpp::NumericVector ncases(VV*GG) ;  
     Rcpp::NumericVector ncases1(VV2) ;  
       
     Rcpp::NumericVector perc1a(VV2) ;  
     Rcpp::NumericVector perc2a(VV2) ;  
       
     Rcpp::NumericMatrix perc1M(VV2,Nimp);  
     Rcpp::NumericMatrix perc2M(VV2,Nimp);  
       
     Rcpp::NumericMatrix perc1_varM(VV2,Nimp);  
     Rcpp::NumericMatrix perc2_varM(VV2,Nimp);  
     Rcpp::NumericMatrix ncases1M(VV2,Nimp);  
       
     Rcpp::NumericMatrix perc1repM(VV2,RR*Nimp);  
     Rcpp::NumericMatrix perc2repM(VV2,RR*Nimp);  
       
       
     Rcpp::Rcout << "|"  ;  
       
     // loop over imputed datasets  
     for (int ii = 0 ;ii<Nimp; ii++ ){   
       
       
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
       
     //****** statistics single data  
     Rcpp::List res1 = bifiehelper_freq( dat1 , wgt1 , group_index1 , group_values ,  
     	vars_values_numb , vars_values , vars_index ,  vars_values_numb_cumsum ) ;  
     perc1a = matr2vec(res1["perc1"]) ;  
     perc2a = matr2vec(res1["perc2"]) ;  
     ncases = res1["ncases"] ;   
     ncases1 = res1["ncases1"] ;  
     sumwgt = matr2vec(res1["sumwgt"]) ;  
       
     //****** statistics replicated dataset  
     Rcpp::List res2 = bifiehelper_freq( dat1 , wgtrep , group_index1 , group_values ,  
     	vars_values_numb , vars_values , vars_index ,  vars_values_numb_cumsum ) ;  
       
     // compute statistics for replicate weights ;  
     Rcpp::NumericMatrix perc1rep = res2["perc1"] ;  
     Rcpp::NumericMatrix perc2rep = res2["perc2"] ;  
       
     // compute standard errors  
     Rcpp::NumericVector perc1_var = varjack_helper( perc1a , perc1rep , fayfac ) ;  
     Rcpp::NumericVector perc2_var = varjack_helper( perc2a , perc2rep , fayfac ) ;	  
       
     for (int zz=0;zz<VV2;zz++){  
          perc1M(zz,ii) = perc1a[zz] ;  
          perc2M(zz,ii) = perc2a[zz] ;  
          perc1_varM(zz,ii) = perc1_var[zz] ;  
          perc2_varM(zz,ii) = perc2_var[zz] ;  
          ncases1M(zz,ii) = ncases1[zz] ;     
          for (int rr=0;rr<RR;rr++){  
     	    perc1repM( zz , rr+ii*RR ) = perc1rep(zz,rr) ;  
     	    perc2repM( zz , rr+ii*RR ) = perc2rep(zz,rr) ;  	      
     	   			   }       
               }  
       
     Rcpp::Rcout << "-" <<  std::flush ;             
       
     } // end ii  
       
     Rcpp::Rcout << "|" << std::endl ;  	  
       
     ///*** Rubin inference  
     Rcpp::List perc1L = rubin_rules_univ( perc1M , perc1_varM ) ;  
     Rcpp::List perc2L = rubin_rules_univ( perc2M , perc2_varM ) ;  
       
     // another output list  
     Rcpp::List out1 = Rcpp::List::create(   
         _["GG"] = GG ,  
         _["VV2"] = VV2  
     		)    ;  
       
       
       
     //*************************************************      
     // OUTPUT              
       
                
     return Rcpp::List::create(   
         _["ncases1M"] = ncases1M ,  
         _["ncases"] = ncases   ,  
         _["perc1"] = perc1L ,  
         _["perc1M"] = perc1M ,  
         _["perc1_varM"] = perc1_varM ,  
         _["perc1repM"] = perc1repM ,  
         _["perc2"] = perc2L ,      
         _["perc2M"] = perc2M   ,  
         _["perc2_varM"] = perc2_varM ,  
         _["perc2repM"] = perc2repM    ,  
         _["outlist"] = out1  
         ) ;    
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;
END_RCPP
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Correlations and covariances


// declarations
extern "C" {
SEXP bifie_correl( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_) ;
}

// definition
SEXP bifie_correl( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
       
     int Nimp = NI[0] ;  
     // int RR = wgtrep.ncol() ;   
       
     int N = wgt1.nrow() ;  
     int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
       
     int WW = wgtrep.ncol() ;  
     int RR = WW ;   
       
     Rcpp::NumericMatrix mean1(VV*GG,WW) ;  
     Rcpp::NumericMatrix sd1(VV*GG,WW) ;  
     Rcpp::NumericVector sumwgt1(GG) ;  
     Rcpp::NumericVector ncases1(GG) ;  
     Rcpp::NumericMatrix ncases1M(GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgt1M(GG,Nimp) ;  
       
       
     // create index matrix for covariances and correlations  
     int ZZ = VV*(VV-1) / 2 + VV ;  
     Rcpp::NumericMatrix itempair_index( ZZ , 2 ) ;  
     Rcpp::NumericVector cov1(ZZ*GG) ;  
     Rcpp::NumericVector cor1(ZZ*GG) ;  
       
     int zz=0;  
     for (int vv1=0;vv1<VV;vv1++){  
       for (int vv2=vv1;vv2<VV;vv2++){	  
          itempair_index(zz,0) = vv1 ;  
          itempair_index(zz,1) = vv2 ;  
          zz++ ;  
        			}  
                }  
       
                  
     int VV2 = ZZ*GG ;             
     Rcpp::NumericMatrix cov1M(VV2,Nimp);  
     Rcpp::NumericMatrix cov1_varM(VV2,Nimp);  
     Rcpp::NumericMatrix cov1repM(VV2,RR*Nimp);  
     Rcpp::NumericMatrix cor1M(VV2,Nimp);  
     Rcpp::NumericMatrix cor1_varM(VV2,Nimp);  
     Rcpp::NumericMatrix cor1repM(VV2,RR*Nimp);  
                  
       
     Rcpp::Rcout << "|"  ;  
       
     ///***************** loop imputed datasets  
                
     for (int ii = 0 ; ii < Nimp ; ii++){   
     // extract dataset  
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
       
     // perform listwise deletion of cases  
     // create indicators of persons  
     Rcpp::NumericVector ind_cases(N) ;   
     for ( int nn=0; nn<N; nn++){ // beg nn  
       ind_cases[nn] = 1 ;  
       for ( int vv=0; vv<VV; vv++){ // beg vv  
           if ( R_IsNA( dat1(nn, vars_index[vv] ) ) ){   
           			ind_cases[nn] = 0 ;   
           				}  
        			}  // end vv  
                  } // end nn	  
       
     // statistics one dataset               
     Rcpp::List res2 = bifiehelpers_correl( dat1 , ind_cases , group_values , group_index1 ,  
     	 wgt1 ,  vars_index , itempair_index ) ;  
     // extract statistics               
     cov1 = matr2vec( res2["cov1"] ) ;               
     cor1 = matr2vec( res2["cor1"] ) ;               
     sumwgt1 = res2["sumwgt1"]  ;               
     ncases1 =  res2["ncases1"]  ;  
       
     // statistics replicated datasets  
     Rcpp::List res3 = bifiehelpers_correl( dat1 , ind_cases , group_values , group_index1 ,  
     	 wgtrep ,  vars_index , itempair_index ) ;  
       
     // compute statistics for replicate weights ;  
     Rcpp::NumericMatrix cov1rep = res3["cov1"] ;  
     Rcpp::NumericMatrix cor1rep = res3["cor1"] ;  
       
     // compute standard errors  
     Rcpp::NumericVector cov1_var = varjack_helper( cov1 , cov1rep , fayfac ) ;  
     Rcpp::NumericVector cor1_var = varjack_helper( cor1 , cor1rep , fayfac ) ;  
       
     for (int zz=0;zz<VV2;zz++){  
          cov1M(zz,ii) = cov1[zz] ;  
          cov1_varM(zz,ii) = cov1_var[zz] ;  
          cor1M(zz,ii) = cor1[zz] ;  
          cor1_varM(zz,ii) = cor1_var[zz] ;          
          for (int rr=0;rr<RR;rr++){  
     	    cov1repM( zz , rr+ii*RR ) = cov1rep(zz,rr) ;  
     	    cor1repM( zz , rr+ii*RR ) = cor1rep(zz,rr) ;  	      
     	   			   }       
               }  
     for (int gg=0;gg<GG;gg++){           
          ncases1M(gg,ii) = ncases1[gg] ;  
          sumwgt1M(gg,ii) = sumwgt1[gg] ;       
          		}  
       
     Rcpp::Rcout << "-" <<  std::flush ;            		  
          		  
          }  // end ii ;  end multiple imputations  
       
     Rcpp::Rcout << "|" << std::endl ;  	       
       
       
     ///*** Rubin inference  
     Rcpp::List cov1L = rubin_rules_univ( cov1M , cov1_varM ) ;  
     Rcpp::List cor1L = rubin_rules_univ( cor1M , cor1_varM ) ;  
       
     // convert output into a matrix (set of matrices)  
     Rcpp::NumericMatrix matr(VV,VV*GG);  
     Rcpp::NumericMatrix matr1(VV,VV*GG);  
       
     //---- estimated correlation  
     Rcpp::NumericVector vec_pars = cor1L["pars"] ;  
     for (int zz=0;zz<ZZ;zz++){  
         int vv1=itempair_index(zz,0) ;	  
         int vv2=itempair_index(zz,1) ;	  
       for (int gg=0;gg<GG;gg++){  
       	 matr( vv1 , vv2 + gg*VV ) = vec_pars[zz*GG + gg ] ;   	    
       	 matr( vv2 , vv1 + gg*VV ) = vec_pars[zz*GG + gg] ;  	    	    	    
            }	  
     }  
     Rcpp::NumericMatrix cor1_matrix = matr ;  
       
       
     //---- estimated correlation  
     vec_pars = cov1L["pars"] ;  
     for (int zz=0;zz<ZZ;zz++){  
         int vv1=itempair_index(zz,0) ;	  
         int vv2=itempair_index(zz,1) ;	  
       for (int gg=0;gg<GG;gg++){  
       	 matr1( vv1 , vv2 + gg*VV ) = vec_pars[zz*GG + gg ] ;   	    
       	 matr1( vv2 , vv1 + gg*VV ) = vec_pars[zz*GG + gg] ;  	    	    	    
            }	  
     }  
     Rcpp::NumericMatrix cov1_matrix = matr1 ;  
          		  
          		  
     //*************************************************      
     // OUTPUT              
       
                
     return Rcpp::List::create(   
         _["itempair_index"] = itempair_index ,  
         _["sumwgt1M"] = sumwgt1M ,  
         _["ncases1M"] = ncases1M ,   
         _["cov1"] = cov1L ,  
         _["cov1M"] = cov1M ,   
         _["cov1repM"] = cov1repM ,  
         _["cov1_varM"] = cov1_varM ,  
         _["cor1"] = cor1L ,  
         _["cor1M"] = cor1M ,   
         _["cor1repM"] = cor1repM ,  
         _["cor1_varM"] = cor1_varM ,  
         _["cor1_matrix"] = cor1_matrix ,  
         _["cov1_matrix"] = cov1_matrix  
         ) ;    
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}



// declarations
extern "C" {
SEXP bifie_linreg( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP dep_index_, 
	SEXP pre_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_) ;
}

// definition

SEXP bifie_linreg( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP dep_index_, 
	SEXP pre_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector dep_index(dep_index_);  
     Rcpp::NumericVector pre_index(pre_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
       
     int Nimp = NI[0] ;  
     // int RR = wgtrep.ncol() ;   
       
     int N = wgt1.nrow() ;  
     int VV = pre_index.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
       
       
     int VV2=(2*VV+2)*GG;  
     Rcpp::NumericMatrix ncasesM(GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgtM(GG,Nimp) ;  
     Rcpp::NumericMatrix regrcoefM(VV2,Nimp) ;  
     Rcpp::NumericMatrix regrcoef_varM(VV2,Nimp) ;  
     int WW = wgtrep.ncol() ;  
     Rcpp::NumericMatrix regrcoefrepM(VV2,Nimp*WW) ;  
       
     Rcpp::Rcout << "|"  ;  
       
     //*****************************  
     // loop over imputations  
     for (int ii = 0 ; ii < Nimp ; ii++){  
       
     // extract dataset  
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
       
     //** apply linear regression to one dataset;  
     Rcpp::List res1 = bifiehelpers_linreg( dat1 , group_values ,  dep_index , pre_index ,   
     	wgt1 , group_index1 )  ;      
     Rcpp::NumericVector ncases = res1["ncases"] ;    
     Rcpp::NumericVector sumwgt0 = matr2vec(res1["sumwgt1"]) ;  
     Rcpp::NumericVector regrcoef0 = matr2vec(res1["regr_coef"]) ;  
       
     //*** apply linear regression to replicated datasets  
     Rcpp::List res2 = bifiehelpers_linreg( dat1 , group_values ,  dep_index , pre_index ,   
     	wgtrep , group_index1 )  ;      
     Rcpp::NumericMatrix sumwgtrep = res2["sumwgt1"] ;  
     Rcpp::NumericMatrix regrcoefrep = res2["regr_coef"] ;  
       
     // compute standard errors  
     Rcpp::NumericVector regrcoef_var = varjack_helper( regrcoef0 , regrcoefrep , fayfac ) ;  
       
     for (int gg=0;gg<GG ; gg++){  
     	ncasesM(gg,ii) = ncases[gg] ;  
     	sumwgtM(gg,ii) = sumwgt0[gg] ;	  
     			}  
     for (int zz=0;zz<VV2; zz++){  
     	regrcoefM(zz,ii) = regrcoef0[zz] ;	  
     	regrcoef_varM(zz,ii) = regrcoef_var[zz] ;  
         for (int ww=0;ww<WW;ww++){  
         	   regrcoefrepM(zz, ww + ii*WW ) = regrcoefrep(zz,ww) ;  
         	   		}  
     }  
       
     Rcpp::Rcout << "-" <<  std::flush ;            		  
          		  
          }  // end ii ;  end multiple imputations  
       
       
       
     Rcpp::Rcout << "|" << std::endl ;  	       
       
     ///*** Rubin inference  
     Rcpp::List regrcoefL = rubin_rules_univ( regrcoefM , regrcoef_varM ) ;       
            
            
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         _["ncasesM"] = ncasesM ,  
         _["sumwgtM"] = sumwgtM ,  
         _["regrcoefrepM"] = regrcoefrepM ,  
         _["regrcoefL"] = regrcoefL ,  
         _["regrcoefM"] = regrcoefM ,      
         _["regrcoef_varM"] = regrcoef_varM   
         ) ;    
       
     // maximal list length is 20!  
              
     
END_RCPP
}

//************************************************+
// BIFIE Wald test


// declarations
extern "C" {
SEXP bifie_waldtest( SEXP parsM_, SEXP parsrepM_, SEXP Cdes_, SEXP rdes_, SEXP Ccols_, SEXP fayfac_) ;
}

// definition

SEXP bifie_waldtest( SEXP parsM_, SEXP parsrepM_, SEXP Cdes_, SEXP rdes_, SEXP Ccols_, SEXP fayfac_ ){
BEGIN_RCPP
   
     //            parsM_ = "matrix" , parsrepM_="matrix" , Cdes_="matrix" ,   
     //			rdes_ = "vector" , Ccols_ = "vector"   
        
        
     Rcpp::NumericMatrix parsM(parsM_);          
     Rcpp::NumericMatrix parsrepM(parsrepM_) ;  
     Rcpp::NumericMatrix Cdes(Cdes_) ;  
     Rcpp::NumericVector rdes(rdes_);  
     Rcpp::NumericVector Ccols(Ccols_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
       
       
     // number of involved variables in the test  
     int VV = Ccols.size() ;  
     int Nimp = parsM.ncol() ;  
     double Nimp2 = Nimp ;  
     Nimp2 = Nimp2 + 1e-10 ;  
     int RR = parsrepM.ncol() / Nimp ;  
     int df = Cdes.nrow() ;  
       
     Rcpp::NumericMatrix chi2M(Nimp,2);  
     // Rcpp::NumericMatrix var_w(VV,VV);  
     arma::mat var_w = arma::zeros(VV,VV);  
     arma::mat var_b = arma::zeros(VV,VV);  
     // Rcpp::NumericMatrix var_b(VV,VV);  
     Rcpp::NumericVector parsM_sel(VV) ;  
       
     // Rcpp::NumericMatrix var_w(VV,VV);  
       
     // arma form of the design matrix Cdes  
     arma::mat ACdes = arma::zeros(df,VV);  
     for (int dd=0;dd<df;dd++){  
         for (int vv=0;vv<VV;vv++){  
         	    ACdes(dd,vv) = Cdes(dd , Ccols[vv] ) ;  
         	    		}  
         	    	}  
        	    	    	    	  
     // arma colvec for design vector  
     arma::colvec Ardes = arma::zeros(df,1) ;  
     for (int dd=0;dd<df;dd++){  
         Ardes(dd,0) =  rdes[dd] ;  
         			}  
         	    	  
     double tmp1=0;  
     double tmp2=0;  
       
     for ( int ii=0 ; ii < Nimp ; ii++){  
     	Rcpp::List res1 = bifiehelpers_waldtest(  VV ,  Ccols , parsM , parsrepM ,  
     		 ii ,  RR , fayfac , ACdes ,  Ardes ) ; 			  
     	  
     	Rcpp::NumericMatrix chi2a=res1["chi2"] ;	  
     	chi2M(ii,0) = chi2a(0,0) ;  
     	tmp1 += chi2M(ii,0) ;	  
     	chi2M(ii,1) = sqrt( chi2M(ii,0) ) ;  
     	tmp2 += chi2M(ii,1) ;		  
     	Rcpp::NumericMatrix var_w1 = res1["var_w"] ;  
     	for (int vv1=0;vv1<VV;vv1++){  
     	   for (int vv2=0;vv2<VV;vv2++){  
     	         var_w(vv1,vv2) += var_w1(vv1,vv2) ;  
     	         			}  
     	         		}  
     		}  
       
     // calculate ARIV	  
     double eps=1e-10;  
     double ariv = tmp1 - Nimp * pow( tmp2 / Nimp , 2.0 ) ;  
     ariv = ariv / ( Nimp - 1 + eps ) * ( 1 + 1 / Nimp2 ) ;  
       
     // calculate D2 statistic  
     double D2 = tmp1 / Nimp2 - (Nimp2+1)/(Nimp2-1+eps) * ariv ;  
     D2 = D2 / ( 1 + ariv ) ;  
     // calculate degrees of freedom  
     double df2 = df ;  
     double nu3 = 1000 ;  
     if ( Nimp > 1 ){  
       nu3 = pow( df2 , - 3 / Nimp2 ) * ( Nimp2 - 1 ) *  
                 pow( 1 + 1 / ( ariv + eps ) , 2 ) ;  
                   }  
       
     double p_D2 = Rf_pf( D2 , df , nu3 , FALSE  , FALSE );                 
       
     // calculate covariance matrices  
     for (int vv1=0;vv1<VV;vv1++){  
        for (int vv2=0;vv2<VV;vv2++){  
              var_w(vv1,vv2) = var_w(vv1,vv2) / Nimp2 ;  
     	         		}  
     	         	}  
       
     // means of all parameters	         	  
     for (int vv=0;vv<VV;vv++){  
        for (int ii=0;ii<Nimp;ii++){  
     	parsM_sel[vv] += parsM( Ccols[vv] , ii ) ;  
     				}  
     	parsM_sel[vv] = parsM_sel[vv] / Nimp2 ;  
     			}  
       
       
     // between matrix  
     // parsM( Ccols[vv1] , ii ) )  
       
     for (int vv1=0;vv1<VV;vv1++){  
     for (int vv2=0;vv2<VV;vv2++){  
       for ( int ii=0; ii<Nimp ; ii++){  
           var_b(vv1,vv2) += ( parsM( Ccols[vv1] , ii ) - parsM_sel[vv1] ) *  
          		( parsM( Ccols[vv2] , ii ) - parsM_sel[vv2] ) ;  
         			}  
           var_b(vv1,vv2) = 1 / ( Nimp2 - 1 ) * var_b( vv1 , vv2 ) ;  
           		}  
               }  
         			  
     arma::mat ariv_D1a = arma::mat( var_b * arma::inv( var_w ) ) ;  
     double ariv_D1 = 0 ;  
     for (int vv=0;vv<VV;vv++){  
     	ariv_D1 += ariv_D1a(vv,vv) ;  
     			}		  
     ariv_D1 = ariv_D1 * ( 1 + 1 / Nimp2 ) / df ;  
     arma::mat var_t1 = arma::mat( (1+ariv_D1) * var_w ) ;  
         			  
     // hypothesis matrix  
     arma::mat var_hyp = arma::mat( ACdes * var_t1 * trans( ACdes) ) ;  
     // compute inverse of variance matrix of hypothesis  
     arma::mat var_hypinv = arma::inv( var_hyp ) ;  
     // parameter vector  
     arma::colvec parm_vec= arma::zeros(VV,1);  
     for (int vv=0;vv<VV;vv++){  
     		parm_vec(vv,0) = parsM_sel[vv] ;  
     				}  
     // hypothesis statistic  
     arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes ) ;  
     arma::mat D1 = arma::mat( trans( hyp_stat ) * var_hypinv * hyp_stat ) ; 			  
     // D1(0,0) = D1(0,0) / df ;  
     // according to Enders (2010, p. 236), D1 must be divided by df,  
     // but this is (could be) an error?  
       
     // calculate nu2  
     double nu2 = 1 + ( 1 - 2 / ( df * Nimp2 - df ) * 1 / ariv_D1 ) ;  
     nu2 = 4 + ( df * Nimp2 - df - 4 ) * nu2 * nu2 ;  

     double tmp11 = D1(0,0);
     
     double p_D1 = Rf_pf( tmp11 , df , nu2 , FALSE  , FALSE );     
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         _["chi2M"] = chi2M ,  
         _["ariv"] = ariv ,  
         _["D2"] = D2 ,  
         _["df"] = df , _["nu2"] = nu2 ,  
         _["nu3"]=nu3 ,  
         _["p_D1"] = p_D1 ,      
         _["p_D2"] = p_D2 ,  
         _["Nimp"] = Nimp , _["RR"] = RR ,  
         _["fayfac"] = fayfac ,  
         _["var_w"] = var_w , _["var_b"] = var_b ,       
         _["D1"] = D1 ,   
         _["Ccols"] = Ccols , _["parsM_sel"] = parsM_sel  
         ) ;    
     // maximal list length is 20!  
              
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}


//*************************************************
// statistical test means


// declarations
extern "C" {
SEXP bifie_test_univar( SEXP mean1M_, SEXP sd1M_, SEXP sumweightM_, SEXP GG_, 
	SEXP group_values_, SEXP mean1repM_, SEXP sd1repM_, SEXP sumweightrepM_, SEXP fayfac_) ;
}

// definition

SEXP bifie_test_univar( SEXP mean1M_, SEXP sd1M_, SEXP sumweightM_, SEXP GG_, 
	SEXP group_values_, SEXP mean1repM_, SEXP sd1repM_, SEXP sumweightrepM_, SEXP fayfac_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix mean1M(mean1M_);  
     Rcpp::NumericMatrix sd1M(sd1M_);  
     Rcpp::NumericMatrix sumweightM(sumweightM_);  
     int GG = as<int>(GG_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
     Rcpp::NumericMatrix mean1repM(mean1repM_);  
     Rcpp::NumericMatrix sd1repM(sd1repM_);  
     Rcpp::NumericMatrix sumweightrepM(sumweightrepM_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
       
       
     int HH = mean1M.nrow() ;  
     int VV = HH / GG ;   
     int GG2 = GG * (GG-1) / 2 ;  
     int Nimp = sd1M.ncol() ;  
       
     Rcpp::NumericMatrix dstatM(VV*GG2,Nimp) ;  
     Rcpp::NumericMatrix dstat_varM(VV*GG2,Nimp) ;  
     Rcpp::NumericMatrix eta2M(VV,Nimp) ;  
     Rcpp::NumericMatrix eta2_varM(VV,Nimp) ;  
     Rcpp::NumericVector eta2V(1) ;  
       
     // matrix of group values  
     Rcpp::NumericMatrix group_values_matrix(GG2,2);  
     int ii=0;  
     for ( int gg1=0; gg1 < GG - 1 ; gg1++){  
     for (int gg2=gg1+1; gg2 < GG ; gg2++){  
     	group_values_matrix(ii,0) = group_values[gg1] ;  
     	group_values_matrix(ii,1) = group_values[gg2] ;	  
     		   ii++ ;  
     			}  
     		}  
       
       
       
     int RR = sd1repM.ncol() / Nimp ;  
       
     Rcpp::Rcout << "|"  ;  
       
     // loop over imputations  
       
     for ( int ii=0; ii < Nimp ; ii++){  
       
     for (int vv=0; vv < VV ; vv++){  
       
       
     // dataset ii  
     Rcpp::NumericMatrix mean1M_ii = mean1M( Range(vv*GG,vv*GG+GG-1) , Range(ii,ii) ) ;  
     Rcpp::NumericMatrix sd1M_ii = sd1M( Range(vv*GG,vv*GG+GG-1) , Range(ii,ii) ) ;  
     Rcpp::NumericMatrix sumweightM_ii = sumweightM( Range(0,GG-1) , Range(ii,ii) ) ;  
       
       
     Rcpp::List res = bifiehelpers_etasquared( mean1M_ii , sd1M_ii , sumweightM_ii ,  GG ) ;  
     eta2V = matr2vec( res["eta2"] ) ;  
     Rcpp::NumericVector dstatV = matr2vec( res["dstat"] ) ;  
       
       
     // Rcpp::NumericMatrix eta2_temp = res["eta2"]  ;  
     // Rcpp::Rcout << "res['eta2']" << eta2_temp(0,0) << std::flush << std::endl ;  
       
       
       
     // analysis replicate weights  
     Rcpp::NumericMatrix mean1M_rr = mean1repM( Range(vv*GG,vv*GG+GG-1) , Range(ii*RR,ii*RR + RR-1) ) ;  
     Rcpp::NumericMatrix sd1M_rr = sd1repM( Range(vv*GG,vv*GG+GG-1) , Range(ii*RR,ii*RR + RR-1)  ) ;  
     Rcpp::NumericMatrix sumweightM_rr = sumweightrepM( Range(0,GG-1) , Range(ii*RR,ii*RR + RR-1)  ) ;  
     Rcpp::List res1 = bifiehelpers_etasquared( mean1M_rr , sd1M_rr , sumweightM_rr ,  GG ) ;  
     Rcpp::NumericMatrix eta2repM = res1["eta2"] ;  
     Rcpp::NumericMatrix dstatrepM = res1["dstat"] ;  
       
     // compute standard errors  
     // eta squared  
     Rcpp::NumericVector eta2_var = varjack_helper( eta2V , eta2repM , fayfac ) ;  
     // d statistics  
     Rcpp::NumericVector dstat_var = varjack_helper( dstatV , dstatrepM , fayfac ) ;  
     // adjusted eta squared statistic  
     // Rcpp::List res3 = varjack_bias_helper( eta2V , eta2repM , fayfac ) ;  
     // Rcpp::NumericVector eta2adj_tmp = res3["pars_bias"] ;  
     // Rcpp::NumericVector eta2adj_var = res3["pars_var"] ;  
     // double eta2adj = eta2V[0] + (RR-1)*( eta2adj_tmp[0] - eta2V[0] ) ;   
     // double eta2adj = RR * eta2V[0] - (RR-1)* eta2adj_tmp[0]  ;  
       
     eta2M(vv,ii) = eta2V[0] ;  
     eta2_varM(vv,ii) = eta2_var[0] ;  
     for (int zz=0;zz<GG2;zz++){  
     	dstatM(zz+vv*GG2,ii) = dstatV[zz] ;  
     	dstat_varM(zz+vv*GG2,ii) = dstat_var[zz] ;  
     			}  
     			  
     			  
     			  
     		} // end vv			  
     Rcpp::Rcout << "-" <<  std::flush ;            		  
          		  
          }  // end ii ;  end multiple imputations  
       
     Rcpp::Rcout << "|" << std::endl ; 			  
     //----------------------  
     		  
     ///*** Rubin inference  
     Rcpp::List eta2L = rubin_rules_univ( eta2M , eta2_varM ) ;       
     Rcpp::List dstatL = rubin_rules_univ( dstatM , dstat_varM ) ;     
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         _["eta2L"] = eta2L ,   
         _["eta2M"] = eta2M ,  
         _["eta2_varM"] = eta2_varM ,  
         _["dstatL"] = dstatL ,       
         _["dstatM"] = dstatM ,  
         _["dstat_varM"] = dstat_varM    ,  
         _["group_values_matrix"] = group_values_matrix  
         ) ;    
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// bifie_crosstab

// declarations
extern "C" {
SEXP bifie_crosstab( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_values1_, 
	SEXP vars_index1_, SEXP vars_values2_, SEXP vars_index2_, SEXP fayfac_, 
	SEXP Nimp_, SEXP group_index_, SEXP group_values_) ;
}

// definition

SEXP bifie_crosstab( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_values1_, 
	SEXP vars_index1_, SEXP vars_values2_, SEXP vars_index2_, SEXP fayfac_, 
	SEXP Nimp_, SEXP group_index_, SEXP group_values_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_values1(vars_values1_) ;  
     Rcpp::NumericVector vars_index1(vars_index1_);  
     Rcpp::NumericVector vars_values2(vars_values2_) ;  
     Rcpp::NumericVector vars_index2(vars_index2_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
       
       
     int Nimp = NI[0] ;  
     int RR = wgtrep.ncol() ;   
     int N = wgt1.nrow() ;  
     int VV1 = vars_values1.size() ;  
     int VV2 = vars_values2.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
     // int WW=1;  
     int ZZ = VV1*VV2*GG ;  
       
     // design matrix  
     Rcpp::NumericMatrix design_pars(ZZ,5) ;  
     int zz=0;  
     for (int gg=0;gg<GG;gg++){  
     for (int vv1=0;vv1<VV1;vv1++){  
     for (int vv2=0;vv2 <VV2;vv2++){  
         design_pars(zz,0) = vars_values1[vv1] ;  
         design_pars(zz,1) = vars_values2[vv2] ;  
         design_pars(zz,2) = group_values[gg] ;  
         design_pars(zz,3) = vv1 ;  
         design_pars(zz,4) = vv2 ;      
         zz++;  
         }  // end vv2  
     }  // end vv1  
     }  // end gg  
       
       
     int CTP = 3*ZZ + VV1*GG + VV2*GG + 2*GG + GG + 3*GG + 3*GG ;  
       
     Rcpp::NumericMatrix ncasesM( ZZ , Nimp ) ;  
     Rcpp::NumericMatrix ncases_ggM( GG , Nimp ) ;  
     Rcpp::NumericMatrix sumwgtM( ZZ , Nimp) ;  
     Rcpp::NumericMatrix ctparsM( CTP , Nimp) ;  
     Rcpp::NumericMatrix ctpars_varM( CTP , Nimp) ;  
     Rcpp::NumericMatrix ctparsrepM( CTP , Nimp*RR) ;  
       
     Rcpp::Rcout << "|"  ;  
       
     ///****** loop imputed datasets  
       
     for (int ii=0; ii <Nimp ; ii++){  
     	  
     	dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
     	  
     	//*** analysis for original data   
     	Rcpp::List res20 = bifiehelpers_crosstab( dat1 ,  wgt1 , group_values , group_index1 ,   
     		vars_values1 ,  vars_index1 , vars_values2 , vars_index2 , design_pars ) ;   
     	Rcpp::NumericMatrix ncases = res20["ncases"] ;  
     	Rcpp::NumericMatrix ncases_gg = res20["ncases_gg"] ;  
     	Rcpp::NumericMatrix sumwgt = res20[ "sumwgt"] ;  
     	Rcpp::NumericMatrix sumwgt_gg = res20["sumwgt_gg"] ;  
     	Rcpp::NumericVector ctpars = matr2vec(res20["crosstab_pars"]);   
     	  
     	  
     	//*** analysis for replicated datasets  
     	Rcpp::List res3 = bifiehelpers_crosstab( dat1 ,  wgtrep , group_values , group_index1 ,   
     		vars_values1 ,  vars_index1 , vars_values2 , vars_index2 , design_pars ) ;   
     	Rcpp::NumericMatrix ctparsrep = res3["crosstab_pars"] ;  
     	  
     	// compute standard errors  
     	Rcpp::NumericVector ctpars_var = varjack_helper( ctpars , ctparsrep , fayfac ) ;  
     	  
     	// data handling  
     	ncasesM(_,ii) = ncases(_,0);  
     	ncases_ggM(_,ii) = ncases_gg(_,0) ;  
     	sumwgtM(_,ii) = sumwgt(_,0);  
     	ctparsM(_,ii) = ctpars ;  
     	ctpars_varM(_,ii) = ctpars_var ;  
     	for (int rr=0;rr<RR;rr++){  
     		ctparsrepM(_,rr+ii*RR) = ctparsrep(_,rr) ;  
     				}  
     	  
     	Rcpp::Rcout << "-" <<  std::flush ;    			  
     		}  
     		  
     Rcpp::Rcout << "|" << std::endl ;  	  
       
     ///*** Rubin inference  
     Rcpp::List ctparsL = rubin_rules_univ( ctparsM , ctpars_varM ) ;  
       
       
     			  
     			  
     //*************************************************      
     // OUTPUT              
                
     return Rcpp::List::create(   
         _["design_pars"] = design_pars ,  
         _["ncases_ggM"] = ncases_ggM ,   
         _["ncasesM"] = ncasesM ,  
         _["sumwgtM"] = sumwgtM ,  
         _["ctparsL"] = ctparsL ,  
         _["ctparsM"] = ctparsM ,  
         _["ctparsrepM"] = ctparsrepM ,  
         _["ctpars_varM"] = ctpars_varM  
         ) ;    
     // maximal list length is 20!  
       
     // Rcpp::Rcout << "(zzi,zzj) " <<  zzi << " " << zzj <<  " vv1=" << vv1 << " vv2=" << vv2 <<   
     //			" h3 =" << h3 << std::flush << std::endl ;  
       
     // Rcpp::Rcout << "(l1x l1y " <<  l1x  << l1y << std::flush << std::endl ;  
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;
END_RCPP
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//* syvby


// declarations
extern "C" {
SEXP bifie_by( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP userfct_) ;
}

// definition

SEXP bifie_by( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP userfct_ ){
BEGIN_RCPP
          
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
       
     Rcpp::Function userfct(userfct_);  
       
       
     int Nimp = NI[0] ;  
     // int RR = wgtrep.ncol() ;   
       
     int N = wgt1.nrow() ;  
     int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
       
     int WW = wgtrep.ncol() ;  
       
       
     // start with a calculation to compute the number of parameters  
     Rcpp::NumericVector w = wgt1(_,0) ;    
     int ii=0;  
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
     Rcpp::NumericMatrix X(N,VV) ;  
     for (int vv=0;vv<VV;vv++){  
         X(_,vv) = dat1(_, vars_index[vv] ) ;	  
     }  
     Rcpp::NumericVector pars= userfct(X,w);  
     int NP = pars.size() ;  
       
       
     Rcpp::NumericMatrix ncasesM(GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgtM(GG,Nimp) ;  
       
     Rcpp::NumericMatrix parsM(NP*GG,Nimp) ;  
     Rcpp::NumericMatrix pars_varM(NP*GG,Nimp) ;  
     Rcpp::NumericVector pars1(NP*GG) ;  
     Rcpp::NumericMatrix pars1rep(NP*GG,WW) ;  
     Rcpp::NumericMatrix pars1repM(NP*GG,WW*Nimp) ;  
     Rcpp::NumericVector pars1_var(NP*GG) ;  
       
     Rcpp::Rcout << "|"  ;  
       
     // dataset ii  
     for (ii=0;ii<Nimp;ii++){  // beg ii  
     	dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;  
     	  
     	//-- compute dimensions  
     	for (int nn=0;nn<N;nn++){ // beg nn  
     	for (int gg =0;gg<GG;gg++){  // beg gg  
     	if ( dat1(nn,group_index) == group_values[gg] ){  
     		ncasesM(gg,ii) ++ ;  
     		sumwgtM(gg,ii) += wgt1(nn,0) ;  
     		break;  
     				}  
     			}  // end gg  
     		} // end nn  
     	  
     	int hh=0;  
     	  
     	for ( int gg=0; gg <GG ; gg++){  // beg gg  
     	   //-- evaluate function	for original dataset ii	  
     		Rcpp::NumericMatrix X1( ncasesM(gg,ii) ,VV) ;  
     		Rcpp::NumericVector w1( ncasesM(gg,ii) ) ;  
     		hh=0;  
     		for (int nn=0;nn<N;nn++){  // beg nn  
     		if ( dat1(nn, group_index ) == group_values[gg] ){ // beg group val  
     			for (int vv=0;vv<VV;vv++){  
     			    X1(hh,vv) = dat1( nn , vars_index[vv] ) ;	  
     				}  
     			w1[hh] = wgt1(nn,0) ;  
     			hh++ ;		  
     			}  // end if group val  
     		}  // end nn  
     		  
     		Rcpp::NumericVector pars_res = userfct( X1 , w1 ) ;  
     		for (int pp=0;pp<NP;pp++){  // beg pp  
     			parsM( pp + gg*NP , ii ) = pars_res[pp] ;  
     			pars1[ pp + gg*NP ] = pars_res[pp] ;	  
     					}  // end pp  
     	    //*** evaluate user function for replicated datasets	  
     		for (int ww=0;ww<WW;ww++){  // beg ww	  
     			Rcpp::NumericVector w2( ncasesM(gg,ii) ) ;  
     			hh=0;  
     			for (int nn=0;nn<N;nn++){  // beg nn  
     			 if ( dat1(nn, group_index ) == group_values[gg] ){  
     				w2[hh] = wgtrep(nn,ww) ;  
     				hh++ ;		  
     				}  
     			}  // end nn  
     			Rcpp::NumericVector pars_res2 = userfct( X1 , w2 ) ;  
     			for (int pp=0;pp<NP;pp++){  // beg pp  
     				pars1repM( pp + gg*NP , ww + ii*WW ) = pars_res2[pp] ;  
     				pars1rep( pp + gg*NP , ww ) = pars_res2[pp] ;	  
     						}  // end pp			  
     				} // end ww  
     		}  // end gg  
     	  
     	// compute standard errors  
     	pars1_var = varjack_helper( pars1 , pars1rep , fayfac ) ;  
     	pars_varM(_,ii) = pars1_var ;  
     Rcpp::Rcout << "-" <<  std::flush ;            		  
          		  
          }  // end ii ;  end multiple imputations  
       
     Rcpp::Rcout << "|" << std::endl ;    
       
     ///*** Rubin inference  
     Rcpp::List parsL = rubin_rules_univ( parsM , pars_varM ) ;  
       
     	  
     	  
     //*************************************************      
     // OUTPUT              
       
                
     return Rcpp::List::create(   
         _["WW"] = WW ,  
         _["N"] = N ,  
         _["NP"] = NP ,      
         _["userfct"] = userfct ,  
         _["sumwgtM"] = sumwgtM ,   
         _["parsrepM"] = pars1repM ,   
         _["parsM"] = parsM ,  
         _["pars_varM"] = pars_varM ,  
         _["ncasesM"] = ncasesM  ,  
         _["parsL"] = parsL  
         ) ;    
     // maximal list length is 20!  
 
END_RCPP
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// histogram


// declarations
extern "C" {
SEXP bifie_hist( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP breaks_) ;
}

// definition

SEXP bifie_hist( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, SEXP breaks_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
     Rcpp::NumericVector breaks(breaks_) ;  
       
     int Nimp = NI[0] ;  
     // int RR = wgtrep.ncol() ;   
       
     int N = wgt1.nrow() ;  
     // int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
     int BB=breaks.size() - 1 ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
       
     Rcpp::NumericMatrix countsM(BB*GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgtM(BB*GG,Nimp) ;  
     Rcpp::NumericMatrix ncasesM(GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgt_ggM(GG,Nimp) ;  
     Rcpp::NumericVector counts(BB*GG) ;  
     Rcpp::NumericVector sumwgt(BB*GG) ;  
     Rcpp::NumericVector sumwgt_gg(GG) ;  
       
     int bb=0;  
       
     Rcpp::Rcout << "|"  ;  
       
     for ( int ii=0; ii < Nimp ; ii++ ){  
       
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;  
       
     for (int nn=0; nn<N;nn++){ // beg nn  
     for (int gg=0; gg<GG;gg++){ // beg gg  
     if ( dat1(nn , group_index ) == group_values[gg] ){  // beg if group val  
        if ( ! R_IsNA( dat1( nn , vars_index[0] ) ) ){  // beg non NA  
     	ncasesM(gg,ii) ++ ;  
     	sumwgt_ggM(gg,ii) += wgt1(nn,0) ;  
     	bb = 0;  
             while (bb < BB ){  // beg while bb  
          	    if ( dat1( nn , vars_index[0] ) < breaks[bb+1] ){ // beg if larger     
          	    	    countsM( bb + gg*BB , ii ) ++ ;  
          	    	    sumwgtM( bb + gg*BB , ii ) += wgt1(nn,0) ;  
          	    	    bb = BB ;  
          	    	    		} // end if larger  
          	     bb ++ ;  
         			} // end while bb  
         	break ;  
        	} // if non NA  
     }  // end if group val  
     } // end if gg  
     } // end nn  
       
     for (int hh=0; hh < BB*GG ; hh++){  
     	counts[hh] += countsM(hh,ii) ;   
     	sumwgt[hh] += sumwgtM(hh,ii) ;  
     				}  
     for (int hh=0; hh < GG ; hh++){  
     	sumwgt_gg[hh] += sumwgt_ggM(hh,ii) ;  
     				}  
     				  
     Rcpp::Rcout << "-" <<  std::flush ;             
       
     } // end ii  
       
     Rcpp::Rcout << "|" << std::endl ;    
       
       
     //*********************************  
     // compute statistics  
       
     for (int hh=0; hh < BB*GG ; hh++){  
     	counts[hh] = counts[hh] / Nimp ;   
     	sumwgt[hh] = sumwgt[hh] / Nimp ;  
     				}  
     for (int hh=0; hh < GG ; hh++){  
     	sumwgt_gg[hh] = sumwgt_gg[hh] / Nimp ;  
     				}  
     				  
     // mid points  
     Rcpp::NumericVector mids(BB-1);  
     for (int bb=0;bb<BB-1;bb++){  
     	mids[bb] = ( breaks[bb] + breaks[bb+1] ) / 2.0 ;  
     			}  
       
     // density  
     Rcpp::NumericVector density_vec(BB*GG) ;  
     Rcpp::NumericVector relfreq(BB*GG) ;  
       
     for (int gg=0;gg<GG; gg++){  
     for (int bb=0;bb<BB;bb++){  
           relfreq[bb+gg*BB] = sumwgt[bb+gg*BB] / sumwgt_gg[gg] ;  
           density_vec[bb+gg*BB] = relfreq[bb+gg*BB] / ( breaks[bb+1] - breaks[bb] ) ;  
           			}  
           		}  
     			  
          		  
     //*************************************************      
     // OUTPUT                                     
     return Rcpp::List::create(   
     	_["BB"] = BB , _["breaks"] = breaks , _["mids"] = mids , 	  
         _["sumwgtM"] = sumwgtM ,  
         _["countsM"] = countsM ,  
         _["ncasesM"] = ncasesM ,  
         _["counts"] = counts ,  
         _["sumwgt"] = sumwgt ,  
         _["sumwgt_ggM"] = sumwgt_ggM ,  
         _["sumwgt_gg"] = sumwgt_gg    ,  
         _["relfreq"] = relfreq ,  
         _["density_vec"] = density_vec  
         ) ;    
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}


//*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
// empirical distribution function


// declarations
extern "C" {
SEXP bifie_ecdf( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, 
	SEXP breaks_, SEXP quanttype_, SEXP maxval_) ;
}

// definition

SEXP bifie_ecdf( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP vars_index_, 
	SEXP fayfac_, SEXP Nimp_, SEXP group_index_, SEXP group_values_, 
	SEXP breaks_, SEXP quanttype_, SEXP maxval_ ){
BEGIN_RCPP
          
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector vars_index(vars_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
     Rcpp::NumericVector breaks(breaks_) ;  
     int quanttype = as<int>(quanttype_) ;  
     int maxval = as<int>(maxval_) ;  
       
       
     int Nimp = NI[0] ;  
     // int RR = wgtrep.ncol() ;   
       
     int N = wgt1.nrow() ;  
     int VV = vars_index.size() ;  
     int NV = datalist.ncol();  
     // int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
     int BB=breaks.size() ;  
       
     int ZZ=VV*GG*BB ;  
     Rcpp::NumericMatrix ecdfM(ZZ,Nimp) ;  
     Rcpp::NumericVector ecdfMtemp(ZZ) ;  
     Rcpp::NumericMatrix ncasesM(VV*GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgtM(VV*GG,Nimp) ;  
     Rcpp::NumericMatrix dat1(N,NV) ;  
       
       
       
     for (int ii=0;ii<Nimp;ii++){ // beg dataset ii  
       
     // int ii=0;  
     	dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;  
     	  
     	ecdfMtemp = bifie_helper_ecdf( dat1 , wgt1 , breaks ,  
     		 group_values , group_index1 ,vars_index , ii ,  
     		 ncasesM ,  sumwgtM , maxval , quanttype ) ;  
       
     	for (int zz=0;zz<ZZ;zz++){  
     		ecdfM(zz,ii) = ecdfMtemp[zz] ;  
     				}  
     // Rcpp::Rcout << "ii= " <<  ii <<  std::flush << std::endl ;  
     	} // end dataset ii  
       
     //********************************	       	  
     // average ecdf  
       
     // Rcpp::Rcout << "before average ecdf" << std::flush << std::endl ;  
     Rcpp::NumericVector ecdf(ZZ);  
     for (int zz=0;zz<ZZ;zz++){  
       for (int ii=0;ii<Nimp;ii++){  
     	ecdf[zz] += ecdfM(zz,ii) ;  
     			}  
     	ecdf[zz] = ecdf[zz] / Nimp ;  
         }  
     // Rcpp::Rcout << "after average ecdf" << std::flush << std::endl ;     	  
       
     //*************************************************      
     // OUTPUT                                     
     return Rcpp::List::create(   
     	_["BB"] = BB ,   
     	_["yval"] = breaks ,  
     	_["sumwgtM"] = sumwgtM ,   
     	_["ncasesM"] = ncasesM ,   
     	_["ecdfM"] = ecdfM  ,  
     	_["ecdf"] = ecdf   
     	) ;  
             
     
END_RCPP
}

//**************************************************************
// logistic regression


// declarations
extern "C" {
SEXP bifie_logistreg( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP dep_index_, 
	SEXP pre_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, 
	SEXP group_values_, SEXP eps_, SEXP maxiter_) ;
}

// definition
SEXP bifie_logistreg( SEXP datalist_, SEXP wgt_, SEXP wgtrep_, SEXP dep_index_, 
	SEXP pre_index_, SEXP fayfac_, SEXP Nimp_, SEXP group_index_, 
	SEXP group_values_, SEXP eps_, SEXP maxiter_ ){
BEGIN_RCPP
   
       
     Rcpp::NumericMatrix datalist(datalist_);          
     Rcpp::NumericMatrix wgt1(wgt_) ;  
     Rcpp::NumericMatrix wgtrep(wgtrep_) ;  
     Rcpp::NumericVector dep_index(dep_index_);  
     Rcpp::NumericVector pre_index(pre_index_);  
     Rcpp::NumericVector fayfac(fayfac_) ;  
     Rcpp::NumericVector NI(Nimp_);  
     Rcpp::NumericVector group_index1(group_index_) ;  
     Rcpp::NumericVector group_values(group_values_) ;  
     double eps=as<double>(eps_) ;  
     int maxiter=as<int>(maxiter_);  
       
     int Nimp = NI[0] ;  
     int RR = wgtrep.ncol() ;  
       
     int N = wgt1.nrow() ;  
     int VV = pre_index.size() ;  
     int NV = datalist.ncol();  
     int group_index = group_index1[0] ;  
     int GG=group_values.size() ;  
       
     Rcpp::NumericMatrix dat1(N,NV) ;  
     Rcpp::NumericVector tempvec(N);  
     int uu=0;  
       
     int VV2=VV*GG;  
     Rcpp::NumericMatrix ncasesM(GG,Nimp) ;  
     Rcpp::NumericMatrix sumwgtM(GG,Nimp) ;  
     Rcpp::NumericMatrix regrcoefM(VV2,Nimp) ;  
     Rcpp::NumericMatrix regrcoef_varM(VV2,Nimp) ;  
     int WW = wgtrep.ncol() ;  
     Rcpp::NumericMatrix regrcoefrepM(VV2,Nimp*WW) ;  
     Rcpp::NumericMatrix tempcoefrepM(VV2,WW) ;  
     Rcpp::NumericVector regrcoef0(VV2) ;  
       
     Rcpp::NumericVector beta0(VV);  
       
     Rcpp::Rcout << "|"  ;  
       
       
     // loop over imputed datasets  
     for ( int ii=0; ii < Nimp ; ii ++ ){  
       
       
     // extract dataset  
     dat1 = datalist( Range( ii*N+0 , ii*N+ (N-1) ) , Range(0,NV-1) ) ;   
       
     // loop over group values gg  
     int ind=1;  
       
       
     for (int gg=0; gg < GG ; gg ++ ){  
     uu=0;  
     for (int nn=0;nn<N;nn++){  
     ind = 1 ;  
     if ( dat1(nn,group_index) == group_values[gg] ){  
       
         if ( R_IsNA( dat1(nn,dep_index[0] ) ) ){  
         	    ind = 0 ;  
         	    	}  // end NA dep  
         for (int vv=0;vv<VV;vv++){  
            if ( R_IsNA( dat1(nn,pre_index[vv] ) ) ){  
         	    ind = 0 ;  
         	    	}  // end NA dep  
         }	    	  
         if ( ind > 0 ){	  
     	ncasesM(gg,ii) ++ ;  
     	sumwgtM(gg,ii) += wgt1[nn] ;  
     	tempvec[uu] = nn ;  
     	uu ++ ;  
     		}  
     			} // end group_val == gg  
          }  
       
     // create datasets for logistic regression  
     int ngg=ncasesM(gg,ii) ;  
     Rcpp::NumericMatrix Xt(ngg,VV);  
     Rcpp::NumericVector yt(ngg);  
     Rcpp::NumericVector wgtt(ngg);  
     Rcpp::NumericMatrix wgtrept(ngg,RR);  
     Rcpp::NumericVector wgttemp(ngg);  
       
     for (int tt=0;tt<ngg;tt++){  
     	yt[tt] = dat1( tempvec[tt] , dep_index[0] ) ;  
     	wgtt[tt] = wgt1[ tempvec[tt] ] ;  
     	for (int rr=0;rr<RR;rr++){  
     		wgtrept(tt,rr) = wgtrep( tempvec[tt],rr) ;  
     				}  
     	for (int vv=0;vv<VV;vv++){  
     		Xt(tt,vv)=dat1( tempvec[tt] , pre_index[vv] ) ;  
     				}  
     			} // end cases tt  
       
     			  
     // logistic regression original dataset  
     Rcpp::List res1 = bifie_estlogistic_helper(  yt ,  
     	Xt , wgtt , beta0 , eps , maxiter ) ;  
       
     Rcpp::NumericVector tempcoef=res1["beta"] ;  
     for (int vv=0;vv<VV;vv++){       
     	regrcoefM(vv+gg*VV,ii) = tempcoef[vv] ;  
     		}  
       
     Rcpp::List res2 ;  
     // replicated datasets  
     // tempcoefrepM  
     for (int rr=0;rr<RR;rr++){		  
     for (int tt=0;tt<ngg;tt++){  
     	wgttemp[tt] = wgtrept(tt,rr) ;  
     			}  
     res2 = bifie_estlogistic_helper(  yt ,	Xt ,   
     	wgttemp , tempcoef , eps , maxiter ) ;  
     Rcpp::NumericVector tempcoef2=res2["beta"] ;  
     for (int vv=0;vv<VV;vv++){       
     	tempcoefrepM(vv+gg*VV,rr) = tempcoef2[vv] ;  
     		}  
     } // end rr  
     } // end gg  
       
     for (int zz=0;zz<VV2;zz++){  
        regrcoef0[zz] = regrcoefM(zz,ii);  
        	}  
        	  
     // compute standard errors  
     Rcpp::NumericVector regrcoef_var = varjack_helper(   
     		regrcoef0 , tempcoefrepM , fayfac ) ;  
       
     for (int zz=0;zz<VV2; zz++){  
     //	regrcoefM(zz,ii) = regrcoef0[zz] ;	  
     	regrcoef_varM(zz,ii) = regrcoef_var[zz] ;  
         for (int ww=0;ww<WW;ww++){  
         	   regrcoefrepM(zz, ww + ii*WW ) = tempcoefrepM(zz,ww) ;  
         	   		}  
     }  
       
       
     Rcpp::Rcout << "-" <<  std::flush ;            		  
          		  
          }  // end ii ;  end multiple imputations  
       
       
       
     Rcpp::Rcout << "|" << std::endl ;  	       
       
       
     ///*** Rubin inference  
     Rcpp::List regrcoefL = rubin_rules_univ( regrcoefM , regrcoef_varM ) ;       
            
       
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         _["ncasesM"] = ncasesM ,  
         _["sumwgtM"] = sumwgtM ,  
         _["regrcoefrepM"] = regrcoefrepM ,  
         _["regrcoefL"] = regrcoefL ,  
         _["regrcoefM"] = regrcoefM ,      
         _["regrcoef_varM"] = regrcoef_varM       
         ) ;    
       
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}


//***********************************************************
// fasttable function


// declarations
extern "C" {
SEXP bifie_fasttable( SEXP datavec_) ;
}

// definition
SEXP bifie_fasttable( SEXP datavec_ ){
BEGIN_RCPP   
       
     Rcpp::NumericMatrix datavec(datavec_);          
       
     int N = datavec.nrow() ;  
       
       
     arma::colvec vals_temp(N) ;  
     int ii=0;  
     for (int nn=0;nn<N;nn++){  
     	if ( ! R_IsNA( datavec(nn,0) ) ){  
     		vals_temp(ii,0) = datavec(nn,0) ;  
     		ii ++ ;  
     					}  
     				}  
       
     int N1 = ii-1 ;				  
     				  
     // arma::mat vec_unique = unique(vals_temp);  
     arma::mat vec_sort = arma::sort( vals_temp( arma::span(0,N1) , arma::span(0,0) ) );  
       
     // create result table  
     Rcpp::NumericMatrix tableM(N1,2);  
     ii = 0 ;  
     tableM(ii,0) = vec_sort(0,0);  
     tableM(ii,1) = 1 ;  
       
     for (int nn=1;nn<N1+1;nn++){  
     	if (vec_sort(nn,0) == tableM(ii,0) ){  
     		tableM(ii,1) ++ ;  
     			} else {  
     		ii ++ ;  
     	        tableM(ii,0) = vec_sort(nn,0) ;  
     	        tableM(ii,1) = 1 ;  
     	        	}  
     	        }  
       
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(  
         _["vec_sort"] = vec_sort ,   
         _["tableM"] = tableM ,  
         _["N_unique"] = ii+1  
         ) ;    
       
     // maximal list length is 20!                
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;       
       
     
END_RCPP
}


//***********************
// table function for characters

// declarations
extern "C" {
SEXP bifie_table1_character( SEXP datavec_) ;
}

// definition

SEXP bifie_table1_character( SEXP datavec_ ){
BEGIN_RCPP
   
       
     Rcpp::CharacterVector datavec(datavec_);                 
     int N = datavec.size() ;  
       
     			  
     Rcpp::CharacterVector uii = Rcpp::unique( datavec ) ; 			  
     Rcpp::IntegerVector indii = Rcpp::match( datavec , uii ) ;  
       
     int Nval = uii.size() ;  
       
     Rcpp::NumericVector tableM(Nval);  
       
     for (int nn=0;nn<N;nn++){  
     	tableM[ indii[nn] - 1 ] ++ ; 	  
     }  
       
       
     //*************************************************      
     // OUTPUT              
               
     return Rcpp::List::create(   
         _["table_names"] = uii ,  
         _["tableM"] = tableM  
         ) ;    

END_RCPP
}





     // maximal list length is 20!                
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;                
     // Rcpp::Rcout << "datavec[nn] " <<  datavec[nn] <<  std::flush << std::endl ;  


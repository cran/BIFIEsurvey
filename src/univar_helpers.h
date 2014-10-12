


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


// user includes

//****************************************
// signum function
Rcpp::NumericVector bifie_sign( double x ){
	Rcpp::NumericVector y1(1);	
	if ( x > 0 ){ y1[0] = 1; }
	if ( x < 0 ){ y1[0] = -1 ; }
	return( wrap(y1) );
		}

//****************************************
// converts a 1-column matrix into a vector
Rcpp::NumericVector matr2vec( Rcpp::NumericMatrix matr1){
	int N1=matr1.nrow();
	Rcpp::NumericVector vect1(N1) ;
	for (int zz=0;zz<N1;zz++){
		vect1[zz] = matr1(zz,0) ;
			}
	return(wrap(vect1))	;
		}


//****************************************************
// matrix entries form parsM into a larger matrix
// starting from row zz in pars_full
Rcpp::List matrix_entry( Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix pars_full1 ,
		int vv ){
	int PP=parsM.nrow();
	int WW=parsM.ncol();	
	int zz=vv;	
	int N1 = pars_full1.nrow() ;
	Rcpp::NumericMatrix pars_full( N1, WW );
	pars_full = pars_full1 ;	
	for (int pp=0;pp<PP;pp++){
	    for (int ww=0;ww<WW;ww++){	
	       pars_full( zz,ww) = parsM(pp,ww) ;
				}
	     zz++ ;
	   }
	Rcpp::NumericVector zz2(1);
	zz2[0] = zz;
	return Rcpp::List::create( 
	    _["pars_full"] = pars_full ,
	    _["zz2"] = zz2 
	    	) ;
	    }

//***************************************************
// mean and SD in case of multiple groups
Rcpp::List univar_helper_multiple_V2group( Rcpp::NumericMatrix dat1, 
     Rcpp::NumericMatrix wgt1,  Rcpp::NumericVector vars_index ,
     Rcpp::NumericVector group_index1 , Rcpp::NumericVector group_values ){	

     int group_index = group_index1[0] ;
     int N = dat1.nrow() ;  
     int VV = vars_index.size() ;  
     int GG=group_values.size() ;
     int WW=wgt1.ncol() ;
     
     Rcpp::NumericMatrix sumwgt1(GG*VV,WW);  
     Rcpp::NumericMatrix ncases1(GG*VV,1);  

     // Rcpp::Rcout << "c200 " << std::flush << std::endl ;
     
     // compute sum of weights  
     for (int nn=0;nn<N;nn++){ // beg nn
        for ( int gg=0;gg<GG;gg++){ // beg gg
          for (int vv=0;vv<VV;vv++){ // beg vv
           if ( ! R_IsNA( dat1(nn,vars_index[vv]) ) ){  // beg if IsNA          	  
            if ( dat1(nn,group_index) == group_values[gg] ){ // beg if dat(,gg)==
              for (int hh=0;hh<WW;hh++){  // beg hh	    
               sumwgt1(gg+vv*GG,hh) += wgt1(nn,hh) ;
               				} // end hh
               ncases1(gg+vv*GG,0) ++ ;  
//               break ;                
               	         } // end if dat(,gg) == group_values[gg]
                         } // end if IsNA
          		} // end vv          		
        	    } // end gg
        	}  // end nn
        	
     int VV1=VV*GG ;   	
     // compute means and standard deviations   		  
     Rcpp::NumericMatrix mean1(VV1,WW);  
     Rcpp::NumericMatrix sd1(VV1,WW);

     Rcpp::NumericMatrix mean1vv(GG,WW);  
     Rcpp::NumericMatrix sd1vv(GG,WW);  
       
     for (int vv=0; vv<VV;vv++){
     	for (int hh=0;hh<WW;hh++){ 
     	   for (int gg=0;gg<GG;gg++){
     		mean1vv(gg,hh) = 0  ;  
     		sd1vv(gg,hh) = 0  ;
     				}
     			}
        for (int nn=0;nn<N;nn++){ // begin nn
       if ( ! R_IsNA( dat1(nn,vars_index[vv]) ) ){  // beg if IsNA            	
        for ( int gg=0;gg<GG;gg++){ // begin gg
            if ( dat1(nn,group_index) == group_values[gg] ){     
              for (int hh=0;hh<WW;hh++){	    
            	mean1vv(gg,hh) += wgt1(nn,hh) * dat1(nn,vars_index[vv]) ;  
                sd1vv(gg,hh) += wgt1(nn,hh) * pow( dat1(nn,vars_index[vv]) , 2.0) ;
                			}
                break ;
                		   } // end if	
                             } // end gg
                         } // end if IsNA                             
     			}  // end nn
     	for (int gg=0;gg<GG;gg++){  // begin gg
        for (int hh=0;hh<WW;hh++){ // beg hh    		
           mean1(vv*GG+gg,hh) = mean1vv(gg,hh) / sumwgt1(gg+vv*GG,hh) ;  
           sd1(vv*GG+gg,hh) = sqrt( ( sd1vv(gg,hh) - sumwgt1(gg+vv*GG,hh)*
           	   	pow(mean1(vv*GG+gg,hh),2) ) /( sumwgt1(gg+vv*GG,hh) - 1 ) ) ;
           				}
           			} // end gg
           
        	} // end vv
        	
        	
	return Rcpp::List::create( 
	    _["sumwgt1"] = sumwgt1 ,
	    _["ncases1"] = ncases1 ,
	    _["mean1"] = mean1 ,
	    _["sd1"] = sd1
	    ) ;          	        	
     }	
//********************************************************************





//*************************************************
// compute standard errors using replicated statistics
Rcpp::NumericVector varjack_helper( Rcpp::NumericVector pars , 
	Rcpp::NumericMatrix pars_jack , Rcpp::NumericVector fayfac ){
	int PP=pars.size();
	// use int here in the subfunction
	int RR=pars_jack.ncol();
	Rcpp::NumericVector pars_var(PP) ;	
	double tmp1=0;	
	for (int pp=0; pp <PP ; pp++){
	   tmp1=0;
	   for (int rr=0;rr<RR;rr++){
	      tmp1 += pow( pars_jack(pp,rr) - pars[pp] ,2.0) ;
				}
	   pars_var[pp] = fayfac[0] * tmp1 ;
			}
    return( wrap( pars_var ) ) ;
    	}

//*************************************************
// compute standard errors and bias correction using replicated statistics
Rcpp::List varjack_bias_helper( Rcpp::NumericVector pars , 
	Rcpp::NumericMatrix pars_jack , Rcpp::NumericVector fayfac ){
	int PP=pars.size();
	// use int here in the subfunction
	int RR=pars_jack.ncol();
	Rcpp::NumericVector pars_var(PP) ;
	Rcpp::NumericVector pars_bias(PP) ;		
	double tmp1=0;	
	for (int pp=0; pp <PP ; pp++){
	   tmp1=0;
	   for (int rr=0;rr<RR;rr++){
	      pars_bias[pp] += pars_jack(pp,rr) ;
				}	   
	   pars_bias[pp] = pars_bias[pp] / RR ;   	   
	   for (int rr=0;rr<RR;rr++){
	      tmp1 += pow( pars_jack(pp,rr) - pars_bias[pp] ,2.0) ;
				}
	   pars_var[pp] = fayfac[0] * tmp1 ;
			}
	return Rcpp::List::create( 
	    _["pars_bias"] = pars_bias ,
	    _["pars_var"] = pars_var      
			    ) ;    
    	}    	
    	
//*************************************************
// Rubin's rules for combining estimates
Rcpp::List rubin_rules_univ( Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix pars_varM){

	int NP = parsM.nrow();
	int Nimp = parsM.ncol();	
	Rcpp::NumericVector pars(NP) ;
	Rcpp::NumericVector pars_varWithin(NP) ;
	Rcpp::NumericVector pars_varBetween(NP) ;
	Rcpp::NumericVector pars_se(NP) ;
	Rcpp::NumericVector pars_fmi(NP) ;
	
	double tmp1=0;
	double tmp2=0;
	double tmp3=0;
	double eps=1e-10;
	double Nimp2 = Nimp ;
	
	for ( int pp = 0; pp < NP ; pp++){
	tmp1=0;
	tmp2=0;  // squared variation , variance between
	tmp3=0; // variance within
	for ( int ii = 0; ii < Nimp ; ii++ ){
	   tmp1 += parsM(pp,ii) ;
	   tmp2 += parsM(pp,ii) * parsM(pp,ii) ;
	   tmp3 += pars_varM(pp,ii) ;
			}
	
	pars[pp] = tmp1 / Nimp ;
	pars_varWithin[pp] = tmp3 / Nimp ;   		
	pars_varBetween[pp] = ( tmp2 - Nimp * pow( pars[pp] , 2.0) ) / (Nimp - 1+eps ) ;
	// ARb 2014-09-10:  added "1.0" instead of "1"
	pars_se[pp] = sqrt( pars_varWithin[pp] + ( 1.0 + 1/Nimp2) * pars_varBetween[pp] ) ;
	pars_fmi[pp] = ( 1.0 + 1/Nimp2) * pars_varBetween[pp] / pow(pars_se[pp] + eps,2.0) ;
	  }

return Rcpp::List::create( 
    _["pars"] = pars ,
    _["pars_se"] = pars_se ,     
    _["pars_varWithin"] = pars_varWithin ,
    _["pars_varBetween"] = pars_varBetween ,
    _["pars_fmi"] = pars_fmi 
		    ) ;  	  	  
	 }

//*******************************************************
// subroutine frequency calculation
Rcpp::List bifiehelper_freq( Rcpp::NumericMatrix dat1 , Rcpp::NumericMatrix wgt ,
	Rcpp::NumericVector group_index1 , Rcpp::NumericVector group_values ,
	Rcpp::NumericVector vars_values_numb , Rcpp::NumericMatrix vars_values ,
	Rcpp::NumericVector vars_index , Rcpp::NumericVector vars_values_numb_cumsum ){

	int WW = wgt.ncol() ; 
	int N=dat1.nrow();
	int group_index = group_index1[0] ;
	int GG=group_values.size() ;
	int VV = vars_index.size() ;	
	int VV2= ( vars_values_numb_cumsum[VV] ) * GG ;

	Rcpp::NumericMatrix perc1(VV2,WW) ;
	Rcpp::NumericMatrix perc2(VV2,WW) ;
	Rcpp::NumericMatrix sumwgt(VV*GG,WW) ;
	Rcpp::NumericVector ncases(VV*GG) ;
	Rcpp::NumericVector ncases1(VV2) ;	
	
        int ind=0;	
	
	//******************************************************
	//********** count frequencies *************************
	for (int nn=0; nn<N ; nn++){ // beg nn
	for (int gg=0; gg<GG;gg++){  // beg gg
	 if ( dat1(nn,group_index) == group_values[gg] ){  // beg if group_values
	  for (int vv=0;vv<VV;vv++){   // beg vv  	 
	  if (! R_IsNA(dat1(nn,vars_index[vv]) ) ){ // beg R_IsNA	  
	   for (int aa=0;aa < vars_values_numb[vv] ; aa++){  // beg aa
	   if ( dat1(nn,vars_index[vv]) == vars_values(aa,vv) ){  // beg if vars_values
		ncases[vv*GG+gg] ++ ;   	   
		ind  = vars_values_numb_cumsum[vv]*GG + aa + gg * vars_values_numb[vv]  ;
		ncases1[ind] ++ ;   
		for (int ww=0;ww<WW;ww++){ // beg ww	
		    perc1(ind,ww) += wgt(nn,ww) ;
		    sumwgt(vv*GG+gg,ww) += wgt(nn,ww) ;
			   }  // end ww
		break ;
		    }  // end if vars_values	            
	      }   // end aa
	  }  // end R_IsNA   
	     } // end vv
	   break ;
		}	// end if = group_values
	  } // end gg			
	} // end nn
	
	//******************************************************
	//********** compute percentages ***********************
	for (int gg=0; gg<GG;gg++){  // beg gg
	  for (int vv=0;vv<VV;vv++){   // beg vv  	   
	   for (int aa=0;aa < vars_values_numb[vv] ; aa++){  // beg aa	   
		ind  = vars_values_numb_cumsum[vv]*GG + aa + gg * vars_values_numb[vv]  ;
		for (int ww=0;ww<WW;ww++){ // beg ww	
		    perc2(ind,ww) = perc1(ind,ww) / sumwgt(vv*GG+gg,ww) ; 
			   }  // end ww            
	      }   // end aa
	     } // end vv
	  } // end gg		
    
	
	return Rcpp::List::create( 
	    _["sumwgt"] = sumwgt ,
	    _["ncases1"] = ncases1 ,
	    _["ncases"] = ncases  ,
	    _["perc1"] = perc1 ,
	    _["perc2"] = perc2   ,
	    _["vars_values_numb_cumsum"] = vars_values_numb_cumsum
	    ) ;  
	}
// Rcpp::Rcout << "c200 " << std::flush << std::endl ;  


//******************************************************
// correlation
Rcpp::List bifiehelpers_correl( Rcpp::NumericMatrix dat1 , Rcpp::NumericVector ind_cases ,
	Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 ,
	Rcpp::NumericMatrix wgt , Rcpp::NumericVector vars_index ,
	Rcpp::NumericMatrix itempair_index	
	){

	int WW = wgt.ncol() ; 
	int N = wgt.nrow() ;
	int VV = vars_index.size() ;
	int group_index = group_index1[0] ;
	int GG=group_values.size() ;
	int ZZ = itempair_index.nrow() ;
	
	Rcpp::NumericMatrix mean1(VV*GG,WW) ;
	Rcpp::NumericMatrix sd1(VV*GG,WW) ;
	Rcpp::NumericMatrix sumwgt1(GG,WW) ;
	Rcpp::NumericVector ncases1(GG) ;
	
	// create index matrix for covariances and correlations
	Rcpp::NumericMatrix cov1(ZZ*GG,WW) ;
	Rcpp::NumericMatrix cor1(ZZ*GG,WW) ;
	
	
	for (int nn=0;nn<N ; nn++){  // beg nn
	if ( ind_cases[nn] == 1 ){  // beg ind_cases == 1     
	for (int gg=0;gg<GG;gg++){	
	   if ( dat1(nn, group_index ) == group_values[gg] ){ // beg group == group_values[gg]
	      ncases1[gg] ++ ;	   
	      for (int ww=0;ww<WW;ww++){  // beg ww	   
		// weights  
		sumwgt1( gg ,ww ) += wgt(nn,ww) ;
		//**** loop variables means and SDs
		for (int vv=0;vv<VV;vv++){
		   int tmpi2 = vv*GG + gg ;		
		   mean1( tmpi2 , ww ) += wgt(nn,ww) * dat1( nn , vars_index[vv] ) ;
		   sd1( tmpi2 , ww ) += wgt(nn,ww) * pow(dat1( nn , vars_index[vv] ),2.0) ;	
					}
		
		// covariances
		for (int zz=0;zz<ZZ;zz++){ // beg zz item pairs
		   int vv1=itempair_index(zz,0) ;
		   int vv2=itempair_index(zz,1) ;
		   if (vv1<vv2){
		     int tmpi1 = zz*GG + gg ;	   		
		     cov1(tmpi1 , ww ) += wgt(nn,ww) * dat1( nn , vars_index[vv1] ) * 
			      dat1( nn , vars_index[vv2] ) ;
				 }
				} // end zz (item pairs)
			     } // end ww
		break ; 	
				} // end group == group_values[gg]			
			} // end gg
			  } // end if ind_cases == 1
			}   // end nn
	
			
	// calculate descriptive statistics
	for (int ww=0;ww<WW;ww++){
	for (int vv=0;vv<VV; vv++){
	   for (int gg=0;gg<GG;gg++){
	       int zz=vv*GG+gg ;	   
	       mean1( zz , ww ) = mean1( zz , ww ) / sumwgt1( gg , ww ) ; 	      
	       sd1(zz, ww) = sqrt( ( sd1(zz,ww) - sumwgt1( gg , ww ) * pow(mean1(zz,ww) ,2.0) )/
			     ( sumwgt1( gg , ww ) - 1 ) ) ;
					}
					}       
					
	for (int pp=0;pp<ZZ;pp++){
	    if (  itempair_index(pp,0) == itempair_index(pp,1) ){
	      int ipp = itempair_index(pp,0) ;	    
	      for (int gg=0;gg<GG;gg++){
		cov1( pp*GG+gg , ww ) = pow( sd1( ipp*GG+gg , ww) , 2.0)  ;           	    
					}
				}
			}
	
			
	// calculate covariance and correlation
	for (int zz=0;zz<ZZ;zz++){
	   int vv1=itempair_index(zz,0) ;
	   int vv2=itempair_index(zz,1) ;
	   if (vv1 < vv2 ){ 	   
	   for (int gg=0;gg<GG;gg++){
	       int tmpi1 = zz*GG + gg ;	 	
		cov1(tmpi1,ww) = ( cov1(tmpi1,ww) - sumwgt1(gg,ww) * mean1( vv1*GG + gg , ww) * 
			mean1( vv2*GG + gg , ww ) ) / ( sumwgt1( gg , ww ) - 1 )  ;
				}
			  }
			}
	
	// calculate correlations    	    	
	for (int zz=0 ; zz < ZZ;zz++){
	  int vv1=itempair_index(zz,0) ;
	  int vv2=itempair_index(zz,1) ;
	  for (int gg=0; gg<GG ;gg++){
	      cor1(zz*GG+gg,ww) = cov1(zz*GG+gg,ww) / sd1( vv1*GG+gg,ww) / sd1( vv2*GG+gg,ww) ;
			}
		}
	}
		
	//**** output		
	return Rcpp::List::create( 
	    _["sumwgt1"] = sumwgt1 , 
	    _["mean1"] = mean1 ,
	    _["sd1"] = sd1 , 
	    _["cov1"] = cov1 ,
	    _["cor1"] = cor1 ,
	    _["ncases1"] = ncases1
	    ) ;  
	}


//***********************************************************
// linear regression
Rcpp::List bifiehelpers_linreg( Rcpp::NumericMatrix dat1 , 
	Rcpp::NumericVector group_values , Rcpp::NumericVector dep_index ,
	Rcpp::NumericVector pre_index , Rcpp::NumericMatrix wgt ,
	Rcpp::NumericVector group_index1 ){
		
int WW = wgt.ncol() ;
int N = wgt.nrow() ;
int VV = pre_index.size() ;
int group_index = group_index1[0] ;
int GG=group_values.size() ;	
int VV2=(2*VV+2)*GG;

Rcpp::NumericMatrix sumwgt1(GG,WW);
Rcpp::NumericVector ncases(GG);
Rcpp::NumericMatrix regr_coef(VV2,WW) ;
Rcpp::NumericMatrix indcases(N,GG);

//**************************
// extract usable cases
for (int nn=0;nn<N;nn++){ // beg nn
   for (int gg=0; gg < GG ; gg++ ){  // beg gg
	if ( dat1(nn,group_index) == group_values[gg] ){ 
		indcases(nn,gg)=1 ;
				}
	if (  R_IsNA( dat1(nn, dep_index[0] ) ) ){ // beg R_IsNA dep var
		indcases(nn,gg)=0;
				} // end R_IsNA dep var
	for (int vv=0;vv<VV;vv++){ // beg vv  independent vars
	   if (  R_IsNA( dat1(nn, pre_index[vv] ) ) ){
		indcases(nn,gg)=0;
				}
			} // end vv  independent vars
        if ( indcases(nn,gg) == 1 ){
        	ncases[gg] ++ ;
        	for (int ww=0;ww<WW; ww++){ 
        	    sumwgt1(gg,ww) += wgt(nn,ww) ;
        	    	}
        	break ;
        		}
		}   // end gg
  }   // end nn	


double sig2=0;
double sig2a=0;

for (int gg=0;gg<GG;gg++){
	int igg=0;
	 
	// create design matrices
	int ngg = ncases[gg] ;
	arma::mat X0=arma::zeros( ngg , VV ) ;
	arma::mat X=arma::zeros( ngg , VV ) ;
	arma::colvec y0=arma::zeros( ngg ) ;
	arma::colvec y=arma::zeros( ngg ) ;
	arma::mat wgt0=arma::zeros( ngg , WW) ;
	
	Rcpp::NumericMatrix M_pre(VV,WW) ;
	Rcpp::NumericMatrix SD_pre(VV,WW) ;
	Rcpp::NumericMatrix M_dep(1,WW) ;
	Rcpp::NumericMatrix SD_dep(1,WW) ;
	
	//***** define input matrices
	for (int nn=0;nn < N ; nn++){  // beg nn
		if ( indcases(nn,gg)==1){  // beg indcases gg 
		    y0(igg,0) = dat1(nn, dep_index[0] ) ;
		    for (int vv=0;vv<VV;vv++){  // beg vv
			   X0(igg,vv) = dat1(nn,pre_index[vv] ) ;
					}  // end vv 	    	   		
		    igg ++ ;
			} // end if indcases gg
		  } // end nn
	
	//**** define used matrices
	double wtmp=0;
	
	for ( int ww=0; ww <WW; ww++){ // beg ww
	igg=0;
	for (int nn=0;nn < N ; nn++){  // beg nn
		if ( indcases(nn,gg)==1){  // beg indcases gg
		    wgt0(igg,ww) = wgt(nn,ww);		
		    wtmp = sqrt( wgt(nn,ww)) ;	
		    y(igg,0) = y0(igg,0)*wtmp ;
		    M_dep(0,ww) += y0(igg,0) * wgt(nn,ww) ;
		    SD_dep(0,ww) += y0(igg,0) * y0(igg,0) * wgt(nn,ww)  ;	    
		    for (int vv=0;vv<VV;vv++){  // beg vv
			   X(igg,vv) = X0(igg, vv )*wtmp ;
			   M_pre(vv,ww) += X0(igg,vv) * wgt(nn,ww) ;
			   SD_pre(vv,ww) += X0(igg,vv) * X0(igg,vv) * wgt(nn,ww) ;	           
					}  // end vv
		    igg ++ ;	    
			} // end if indcases gg
		  } // end nn	  
	
	//*** fit linear model	  
	arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
	// arma::colvec resid = y - X*coef;            // residuals
	// not that the weights are already included in the residual calculation
	arma::colvec resid(ngg);
	for (int hh=0;hh<ngg;hh++){
		resid(hh,0) = y0(hh,0) ;
		for (int vv=0;vv<VV;vv++){
		   resid(hh,0) = resid(hh,0) - X0(hh,vv) * coef(vv,0) ;
			}	
	}
	     
	// sig2 = arma::as_scalar( arma::trans(resid)*resid );
	sig2=0;
	for (int hh=0;hh<ngg;hh++){
	      sig2 += pow( resid(hh,0),2.0) * wgt0(hh,ww) ;	
	   }
	double sggww = sumwgt1(gg,ww) ;
	sig2a = sig2 / ( sggww - VV ) ;
	
	// collect all regression coefficients
	// unstandardized coefficients
	for (int vv=0;vv<VV;vv++){
		regr_coef( vv + gg*(2*VV+2) , ww ) = coef(vv,0) ;
			}
	regr_coef( VV + gg*(2*VV+2) , ww ) = sqrt( sig2a ) ;  // sigma		  
	// compute R^2
	M_dep(0,ww) = M_dep(0,ww) / sggww ;
	double sig3 = SD_dep(0,ww) - sggww * pow( M_dep(0,ww) , 2.0) ; 
	SD_dep(0,ww) = sqrt( sig3 / ( sggww - 1 ) ) ;
	regr_coef( VV+1 + gg*(2*VV+2) , ww ) = 1 - sig2 / sig3 ;
	
	// compute standardized coefficients
	for (int vv=0;vv<VV;vv++){
		M_pre(vv,ww) = M_pre(vv,ww) / sggww ;
		SD_pre(vv,ww) = SD_pre(vv,ww) - sggww *  pow( M_pre(vv,ww) , 2.0) ;
		SD_pre(vv,ww) = sqrt( SD_pre(vv,ww) / ( sggww - 1 ) ) ;
		regr_coef( VV+2+vv + gg*(2*VV+2) , ww ) = 
		   coef(vv,0) / SD_dep(0,ww) * SD_pre(vv,ww) ;
				}
	} // end ww
} // end gg
       
return Rcpp::List::create( 
    _["ncases"] = ncases ,
    _["sumwgt1"] = sumwgt1 ,
    _["regr_coef"] = regr_coef 
    ) ; 


}

//********************************************
// BIFIE helpers Wald test
Rcpp::List bifiehelpers_waldtest( int VV , Rcpp::NumericVector Ccols ,
	Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix parsrepM ,
	int ii , int RR , Rcpp::NumericVector fayfac ,
	arma::mat ACdes , arma::colvec Ardes ){
	
	//*** calculate covariance matrix for imputation ii
	arma::mat var_w = arma::zeros(VV,VV);
	for (int vv1=0;vv1<VV;vv1++){
	for (int vv2=vv1;vv2<VV;vv2++){
	for (int rr=0;rr<RR;rr++){
	   var_w(vv1,vv2) += ( parsrepM( Ccols[vv1] , rr+ii*RR ) - parsM( Ccols[vv1] , ii ) )*
			( parsrepM( Ccols[vv2] , rr+ii*RR ) - parsM( Ccols[vv2] , ii ) ) ;
				}	
	var_w(vv1,vv2) = fayfac[0] * var_w(vv1,vv2) ;
	var_w(vv2,vv1) = var_w(vv1,vv2) ;
				}
			}
	
	//***  compute chi squared Wald test statistic
	// compute covariance matrix of hypothesis
	arma::mat var_hyp = arma::mat( ACdes * var_w * trans( ACdes) ) ;
	
	//@@@@@ ARb 2014-07-31
//	int XX= arma::mat( var_hyp.n_rows ) ;
//	double eps2 = 1e-15;
//	for (int xx=0;xx<XX;xx++){
//		var_hyp(xx,xx) = var_hyp(xx,xx) + eps2 ;
//				}	
	//@@@@@
	
	// compute inverse of variance matrix of hypothesis
	arma::mat var_hypinv = arma::inv( var_hyp ) ;
	// parameter vector
	arma::colvec parm_vec= arma::zeros(VV,1);
	for (int vv=0;vv<VV;vv++){
		parm_vec(vv,0) = parsM( Ccols[vv] , ii ) ;
				}
	// hypothesis statistic
	arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes ) ;
	arma::mat chi2 = arma::mat( trans( hyp_stat ) * var_hypinv * hyp_stat ) ; 			
	
	return Rcpp::List::create( 
	    _["chi2"] = chi2 ,
	    _["var_w"] = var_w
	    ) ; 
	}

	
//*******************************************************
// etasquared and d statistics
Rcpp::List bifiehelpers_etasquared( Rcpp::NumericMatrix mean1M ,
	Rcpp::NumericMatrix sd1M , Rcpp::NumericMatrix sumweightM , int GG ){
	
	int WW = sd1M.ncol() ;
	// calculate total mean
	Rcpp::NumericVector totmean(WW) ; 
	Rcpp::NumericVector sumwgt(WW) ; 
	Rcpp::NumericVector expl_var(WW) ;
	Rcpp::NumericVector resid_var(WW) ;
	Rcpp::NumericMatrix eta2(1,WW) ;
	int GG2 = GG * (GG - 1 ) / 2 ;
	Rcpp::NumericMatrix dstat(GG2 , WW ) ;
	

	for ( int ww=0; ww < WW ; ww++){ // beg ww
		for (int gg=0;gg<GG; gg++){ // beg gg
		  sumwgt[ww] += sumweightM(gg,ww) ;
		  totmean[ww] += sumweightM(gg,ww) * mean1M(gg,ww) ;	  
// Rcpp::Rcout << "mean1M(gg,0)=" << mean1M(gg,ww) << std::flush << std::endl ;

			}  // end gg
		totmean[ww] = totmean[ww] / sumwgt[ww] ;  
		for (int gg=0;gg<GG; gg++){  // beg gg
		  expl_var[ww] += sumweightM(gg,ww)*pow( mean1M(gg,ww) - totmean[ww] , 2.0 ) ;
		  resid_var[ww] += (sumweightM(gg,ww)-1)*pow( sd1M(gg,ww) , 2.0 ) ; 
		  eta2(0,ww) = sqrt( expl_var[ww] / ( expl_var[ww] + resid_var[ww] ) ) ;
			}  // end gg
		
		// calculate d statistics
		int ii=0;
		for ( int gg1=0; gg1 < GG - 1 ; gg1++){
		for (int gg2=gg1+1; gg2 < GG ; gg2++){
		   dstat(ii,ww) = mean1M(gg1,ww) - mean1M(gg2,ww) ;
		   dstat(ii,ww) = dstat(ii,ww) / sqrt( 0.5 * ( pow(sd1M(gg1,ww),2.0) + 
				    pow(sd1M(gg2,ww),2.0) ) ) ;
		   ii++ ;
				}
			}
	 }  // end ww


	return Rcpp::List::create( 
	    _["eta2"] = eta2 ,
	    _["dstat"] = dstat
	    ) ;   
  }
  
  
//***********************************************************
// cross tabulation
Rcpp::List bifiehelpers_crosstab( Rcpp::NumericMatrix dat1 , Rcpp::NumericMatrix wgt ,
	Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 , 
	Rcpp::NumericVector vars_values1 , Rcpp::NumericVector vars_index1 ,
	Rcpp::NumericVector vars_values2 , Rcpp::NumericVector vars_index2 ,	
	Rcpp::NumericMatrix design_pars 
		){
   int N=dat1.nrow();   
   int group_index = group_index1[0] ;        
   int GG=group_values.size() ;
   int WW=wgt.ncol();
   int VV1 = vars_values1.size() ;
   int VV2 = vars_values2.size() ;
   int ZZ = VV1*VV2*GG ;
   int zz=0;
   
   Rcpp::NumericMatrix ncases( ZZ , 1 ) ;
   Rcpp::NumericMatrix sumwgt( ZZ , WW ) ;
	for (int nn=0;nn<N; nn++){  // beg nn
	for (int gg=0;gg<GG;gg++){ // beg gg
	if ( dat1( nn , group_index ) == group_values[gg] ){  // beg if groupval gg
	for (int vv1=0;vv1<VV1;vv1++){   // beg vv1
	   for (int vv2=0;vv2<VV2;vv2++){	// beg vv2
	      if ( ( dat1(nn, vars_index1[0]) == vars_values1[vv1] ) &
		       ( dat1(nn, vars_index2[0]) == vars_values2[vv2] ) ){  // beg if comb
			ncases( vv2 + vv1*VV2 + gg*VV1*VV2 , 0 ) ++ ;
			for (int ww=0;ww<WW;ww++){  // beg ww
			   sumwgt( vv2 + vv1*VV2 + gg*VV1*VV2 , ww ) += wgt(nn,ww) ; 
			      } // end ww
                        break ; 			      
				}  // end if combination of values
			}  // end vv2
		}   // end vv1
	    } // end if groupval gg
	} // end gg
	}  // end nn
	
	Rcpp::NumericMatrix ncases_gg(GG,1);
	Rcpp::NumericMatrix sumwgt_gg(GG,WW);
	
	Rcpp::NumericMatrix probs_joint(ZZ,WW) ;
	Rcpp::NumericMatrix probs_rowcond(ZZ,WW) ; 
	Rcpp::NumericMatrix probs_colcond(ZZ,WW) ;
	
	Rcpp::NumericMatrix probs_rowmarg(GG*VV1,WW);
	Rcpp::NumericMatrix probs_colmarg(GG*VV2,WW);
	
	//**** compute counts within each group value gg
	for (int gg=0;gg<GG;gg++){ // beg gg
	for (int vv1=0;vv1<VV1;vv1++){  // beg vv1
	   for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
		   ncases_gg(gg,0) += ncases( vv2 + vv1*VV2 + gg*VV1*VV2 , 0 ) ;
		  for (int ww=0;ww<WW;ww++){
			sumwgt_gg(gg,ww) += sumwgt( vv2 + vv1*VV2 + gg*VV1*VV2 , ww ) ;
					}
					}  // end vv2
				}  // end vv1
		   } // end gg
	
	//**** compute joint probabilities
	int ind=0;
	for (int gg=0;gg<GG;gg++){ // beg gg
	for (int vv1=0;vv1<VV1;vv1++){  // beg vv1
	   for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
		  ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
		  for (int ww=0;ww<WW;ww++){
		     probs_joint(ind,ww)= sumwgt( ind , ww ) / sumwgt_gg( gg ,ww);
					}
				}  // end vv2
		   }  // end vv1
		} // end gg	   
	//**** compute conditional rowwise probabilities   	   
	double tmp1=0;
	for (int ww=0;ww<WW;ww++){
	for ( int gg=0;gg<GG;gg++){
	for (int vv1=0;vv1<VV1;vv1++){ // beg vv1
		tmp1=0;
		for (int vv2=0;vv2<VV2;vv2++){
		   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
		   tmp1 += sumwgt( ind , ww ) ;
					}
		   probs_rowmarg( vv1 + gg*VV1 , ww ) = tmp1 / sumwgt_gg( gg , ww ) ;  
		for (int vv2=0;vv2<VV2;vv2++){
		   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
		   probs_rowcond( ind , ww ) = sumwgt( ind , ww ) / tmp1 ;
					}
	    } // end vv1	
	} //end gg
	} // end ww
	
	//**** compute conditional columnwise probabilities
	for (int ww=0;ww<WW;ww++){
	for (int gg=0;gg<GG;gg++){
	for (int vv2=0;vv2<VV2;vv2++){
	tmp1=0;
	     for (int vv1=0;vv1<VV1;vv1++){
		ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
		tmp1 += sumwgt( ind , ww ) ;
				} // end vv1
	     probs_colmarg( vv2 + gg*VV2 , ww ) = tmp1 / sumwgt_gg(gg,ww) ;
	     for (int vv1=0;vv1<VV1;vv1++){
		   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
		   probs_colcond( ind , ww ) = sumwgt( ind , ww ) / tmp1 ;
					}			
	}   // end vv2
	} // end gg
	} // end ww
		   
	//*** effect size w  (also include Cramer's V)
	Rcpp::NumericMatrix w_es(2*GG,WW);  // w and Cramer's V
	int VV3 = VV1 ;
	if (VV2 < VV1 ){ VV3 = VV2 ; }
	double tmp2=0;
	for (int ww=0; ww<WW;ww++){ // beg ww
	for (int gg=0;gg<GG;gg++){ // beg gg
	for (int vv1=0;vv1<VV1;vv1++){ // beg vv1
	for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
	   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
	   // probability under independence;
	   tmp2 = probs_rowmarg( vv1 + gg*VV1 , ww ) * probs_colmarg( vv2 + gg*VV2, ww );
	   w_es(gg,ww) += pow( probs_joint( ind , ww ) - tmp2 , 2.0 ) / tmp2 ;
	//  Rcpp::Rcout << "w_es " <<  w_es(gg,0) << "  gg,ww" << gg << "  " << ww <<  std::flush << std::endl ;   
					} // end vv2
	     } // end vv1
	w_es(gg,ww) = sqrt( w_es(gg,ww) ) ;  // effect size w
	w_es(GG+gg,ww ) = w_es(gg,ww) / sqrt( VV3 - 1.0 ) ; // Cramer's V
	} // end gg
	} // end ww
	
	//*** gamma measure of association (Goodman)
	Rcpp::NumericMatrix gamma_es( GG , WW ) ;
	
	double pc=0;
	double pd=0;
	// double h1=0;
	double h2=0;
	int VV12 = VV1*VV2 ;
	
	for (int ww=0;ww<WW;ww++){   // beg ww
	for (int gg=0; gg < GG ; gg++){ // beg gg
		pc=0;
		pd=0;
		// numer := 0
		//for i:=2..N do
		//    for j:=1..(i-1) do
		//        numer := numer + sign(x[i] - x[j]) * sign(y[i] - y[j])
		// return numer
		
		for (int zzi = 0 ; zzi < VV12 ; zzi++ ){ // beg zzi
		for( int zzj = 0 ; zzj < VV12 ; zzj++ ){ // beg zzj
		
		double h1 = vars_values1[ design_pars(zzi , 3) ] - vars_values1[ design_pars(zzj , 3) ] ;
		Rcpp::NumericVector h1s = bifie_sign( h1 ) ;
		
		h2 = vars_values2[ design_pars(zzi , 4) ] > vars_values2[ design_pars(zzj , 4) ] ; 
		Rcpp::NumericVector h2s = bifie_sign( h2 ) ;
		
		double h3 = h1s[0] * h2s[0] ;
		if ( zzi== zzj ){ h3 = 0 ; }
		int vv1 = design_pars(zzj , 3) ;
		int vv2 = design_pars(zzj , 4) ;
		ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;   
		
		double p3= probs_joint( design_pars(zzi,4) + design_pars(zzi , 3)*VV2 + gg*VV1*VV2 , ww) ;		
				
		// concordant cases
		if ( h3 > .99 ){
		    pc += probs_joint(ind,ww) * p3;
		   }  // end if concordant
		
		// discordant cases
		if ( h3 < -.99 ){ 
		    pd += probs_joint(ind,ww) * p3 ;
		   }  // end if discordant   
		}  // end zzj
		} // end zzi
		gamma_es(gg,ww) = ( pc - pd ) / ( pc + pd ) ;
	} // end gg
	} // end ww
	//**** lambda PRE measures
	Rcpp::NumericMatrix lambda(3*GG, WW) ;  // lambda_V1, lambda_V2, lambda
	double l1y=0;
	double l2y=0;
	double p1=0;
	double t11=0;
	double l1x=0;
	double l2x=0;
	
	for (int ww=0; ww<WW;ww++){ // beg ww
	for (int gg=0;gg<GG;gg++){ // beg gg
	
		//-- lambda_Y
		l1y=0;
		l2y=0;
		for (int vv1 = 0 ;vv1<VV1;vv1++){
			t11 = 0 ;
			for (int vv2=0;vv2<VV2;vv2++){
			     ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
			     p1 = probs_joint( ind , ww ) ;	
			     if (p1> t11 ){ t11 = p1 ; } 
				}
			l1y += t11 ;
			}	
		for (int vv2=0;vv2<VV2;vv2++){
		   p1 = probs_colmarg( vv2 + gg*VV2 , ww ) ;
		   if (p1 > l2y ){ l2y = p1 ; }
					}
		lambda( 2 + gg*3 , ww ) = ( l1y - l2y ) / ( 1 - l2y ) ; 
		
		//-- lambda_X
		l1x=0;
		l2x=0;
		for (int vv2 = 0 ;vv2<VV2;vv2++){
			t11 = 0 ;
			for (int vv1=0;vv1<VV1;vv1++){
			     ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
			     p1 = probs_joint( ind , ww ) ;	
			     if (p1> t11 ){ t11 = p1 ; } 
				}
			l1x += t11 ;
			}	
		for (int vv1=0;vv1<VV1;vv1++){
		   p1 = probs_rowmarg( vv1 + gg*VV1 , ww ) ;
		   if (p1 > l2x ){ l2x = p1 ; }				
		   } 
		lambda( 1 + gg*3 , ww ) = ( l1x - l2x ) / ( 1 - l2x ) ; 
		
		//-- symmetric lambda measure
		lambda( 0 + gg*3 , ww ) = ( l1x + l1y - l2x - l2y ) / ( 2 - l2x - l2y ) ;
	
		
	} // end gg
	} // end ww
	
	//**** Kruskal's tau	
	Rcpp::NumericMatrix kruskal_tau(3*GG,WW) ;
	
	double t1y=0 ;
	double t2y=0 ;
	double t1x=0;
	double t2x=0;
	
	
	for (int ww=0;ww<WW;ww++){ // beg ww
	for (int gg=0;gg<GG;gg++){ // beg gg
	//*** tau_Y
	t1y=0;
	t2y=0;
	for (int vv2=0;vv2<VV2;vv2++){  // beg vv1
		t2y += pow( probs_colmarg( vv2 + gg*VV2 , ww ) , 2.0 );
				  }  // end vv2
	for ( int vv1=0; vv1 < VV1 ; vv1++){  // beg vv1
	for (int vv2=0;vv2<VV2;vv2++){   // beg vv2
	   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
	   t1y += pow( probs_joint( ind , ww ) , 2.0 ) / probs_rowmarg( vv1+gg*VV1, ww ) ;
				}   // end vv2
			}  // end vv1
	kruskal_tau(2+gg*3,ww) = ( t1y - t2y ) / ( 1 - t2y ) ;
	
	//*** tau_X
	t1x=0;
	t2x=0;
	for (int vv1=0;vv1<VV1;vv1++){
		t2x += pow( probs_rowmarg( vv1 + gg*VV1 , ww ) , 2.0 );
				  }
	for ( int vv2=0; vv2 < VV2 ; vv2++){  // beg vv2
	for (int vv1=0;vv1<VV1;vv1++){   // beg vv1
	   ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
	   t1x += pow( probs_joint( ind , ww ) , 2.0 ) / probs_colmarg( vv2+gg*VV2, ww ) ;
				}   // end vv1
			}  // end vv2
	kruskal_tau(1+gg*3,ww) = ( t1x - t2x ) / ( 1 - t2x ) ;
	
	//*** symmetric tau
	kruskal_tau(0+gg*3,ww) = ( t1x+t1y-t2x-t2y)/(2-t2x-t2y);
	}  // end gg
	} // end ww
	
	//**** output vector of parameters from crosstable
	
	// probs_joint    ZZ
	// probs_rowcond  ZZ
	// probs_colcond   ZZ
	// probs_rowmarg   VV1*GG
	// probs_colmarg  VV2*GG
	// w_es  2*GG
	// gamma_es  GG
	// lambda   3*GG
	// kruskal_tau  3*GG
	
	int CTP = 3*ZZ + VV1*GG + VV2*GG + 2*GG + GG + 3*GG + 3*GG ;
	Rcpp::NumericMatrix crosstab_pars( CTP , WW) ;
	zz=0;
	//*** probs_joint
	Rcpp::List res2= matrix_entry( probs_joint , crosstab_pars ,zz ) ;   
	Rcpp::NumericVector zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v2 = res2["pars_full"] ;
	//*** rowcond
	res2= matrix_entry( probs_rowcond , v2 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v3 = res2["pars_full"] ;
	//*** colcond
	res2= matrix_entry( probs_colcond , v3 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v4 = res2["pars_full"] ;
	//*** rowmarg
	res2= matrix_entry( probs_rowmarg , v4 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v5 = res2["pars_full"] ;
	//*** colmarg
	res2= matrix_entry( probs_colmarg , v5 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v6 = res2["pars_full"] ;
	//*** w_es
	res2= matrix_entry( w_es , v6 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v7 = res2["pars_full"] ;
	//*** gamma_es
	res2= matrix_entry( gamma_es , v7 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v8 = res2["pars_full"] ;
	//*** lambda
	res2= matrix_entry( lambda , v8 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v9 = res2["pars_full"] ;
	//*** kruskal_tau
	res2= matrix_entry( kruskal_tau , v9 ,zz ) ;
	zz2 = res2["zz2"] ;
	zz=zz2[0] ;
	Rcpp::NumericMatrix v10 = res2["pars_full"] ;	
	//*************************************************    
	// OUTPUT                    
	return Rcpp::List::create( 
	    _["ncases"] = ncases ,
	    _["ncases_gg"] = ncases_gg , 
	    _["sumwgt"] = sumwgt ,
	    _["sumwgt_gg"] = sumwgt_gg ,
	    _["crosstab_pars"] = v10 
	    ) ;  
    }
    
/////////////////////////////////////////////////////    
// bifie_helper_ecdf
Rcpp::NumericVector bifie_helper_ecdf( Rcpp::NumericMatrix dat1 , 
	Rcpp::NumericVector wgt , Rcpp::NumericVector breaks ,
	Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 ,
	Rcpp::NumericVector vars_index , int ii ,
	Rcpp::NumericMatrix ncasesM , Rcpp::NumericMatrix sumwgtM ,
	int maxval , int quanttype ){
		
	int N=dat1.nrow();	
	int GG=group_values.size() ;
	int BB=breaks.size() ;
	int VV=vars_index.size();
	Rcpp::NumericVector ecdf_temp(BB);
	arma::colvec vals_temp(N) ;
	arma::colvec wgt_temp(N) ;
	arma::colvec group_temp(N) ;
	arma::uvec indvv(N);
	arma::colvec wgt_tempsort(N);
	arma::colvec vals_tempsort(N);
	arma::colvec group_tempsort(N);
	arma::colvec wgt_tempsort_gg(N);
	arma::colvec vals_tempsort_gg(N);
	int group_index = group_index1[0] ;
	int ZZ=VV*BB*GG;
	Rcpp::NumericVector ecdfMtemp(ZZ);
	
	Rcpp::NumericMatrix a1(N,4) ;
	double eps = 1e-20 ;
	double maxval2 = maxval ;
	
	int uu=0;
	int ngg=0;
	int hh=0;
	int MM=0;
	
	//*** count observed responses
	for (int nn=0;nn<N;nn++){ // beg nn
	for (int gg=0;gg<GG;gg++){ // beg gg
	if ( dat1(nn,group_index ) == group_values[gg] ){ // beg val gg
	for (int vv=0;vv<VV;vv++){ // beg vv
	    if ( ! R_IsNA( dat1(nn,vars_index[vv] ) ) ){	
		ncasesM(gg + vv*GG , ii ) ++ ;
		sumwgtM(gg+vv*GG,ii) += wgt(nn,0) ;
					}
				}  // end vv
		}  // enf if val gg
	} // end gg
	} // end nn


	int bb=0;
	
	for (int vv=0;vv<VV;vv++){
	
	for (int nn=0;nn<N;nn++){ // beg nn
	   group_temp(nn,0) = -1 ;
	   wgt_temp(nn,0) = 0 ;
	   vals_temp(nn,0) = maxval2 ;
	   for (int gg=0;gg<GG;gg++){
	   if ( dat1(nn,group_index ) == group_values[gg] ){ // beg val gg
	      if ( ! R_IsNA( dat1(nn,vars_index[vv] ) ) ){ // beg non NA
			   vals_temp(nn,0) = dat1(nn,vars_index[vv] ) ;
			   wgt_temp(nn,0) = wgt(nn,0 ) ; 
			   group_temp(nn,0) = gg ;
				} // end non NA
	     } // end val gg
	   } // end gg
	} // end nn 
		
		
	
	// sort values
	indvv = arma::sort_index( vals_temp ) ;
	vals_tempsort =vals_temp.rows( indvv ) ;
	wgt_tempsort =wgt_temp.rows( indvv ) ;
	group_tempsort =group_temp.rows( indvv ) ;
	
	
	for ( int gg=0;gg<GG;gg++){
	ngg = ncasesM( gg + vv*GG , ii ) ;
	uu=0;
	for (int nn=0;nn<N;nn++){ // beg nn
	   if ( group_tempsort(nn,0) == gg ){	// beg group = gg
	     vals_tempsort_gg(uu,0) = vals_tempsort(nn,0) ;
	     wgt_tempsort_gg(uu,0) = wgt_tempsort(nn,0) / sumwgtM( gg+vv*GG, ii) ;
	     uu ++ ;
	      } // end if group = gg
	}  // end nn
	
	
	if (ngg<N){    	// beg if ngg < N
	    for (int uu=ngg; uu<N;uu++){  // beg uu
		 vals_tempsort_gg(uu,0) = 0 ;
		 wgt_tempsort_gg(uu,0) = 0 ;
			} // end uu
		} // end if ngg < N
		
	// define ECDF	
	a1(0,0)=0;
	a1(0,1)=wgt_tempsort_gg(0,0);
	a1(0,2)=vals_tempsort_gg(0,0);
	a1(0,3)=vals_tempsort_gg(0,0);
	
	for (int uuu=1;uuu<ngg;uuu++){
		a1(uuu,0) = a1(uuu-1,1) ;
		a1(uuu,1) = a1(uuu-1,1) + wgt_tempsort_gg(uuu,0) ;
		a1(uuu,2) = a1(uuu-1,3);
		a1(uuu,3) = vals_tempsort_gg(uuu,0) ;
			}
	
			
	uu=0;		
	bb=0;
	MM=0;
	while ( ( bb < BB ) & ( ! ( uu >= ngg ) ) ){
		if ( a1(uu,1) > breaks[bb] ){
		   if (quanttype==1){	
			ecdf_temp[bb] = a1(uu,2) + ( a1(uu,3) - a1(uu,2) ) * 
						( breaks[bb] - a1(uu,0)  ) / ( a1(uu,1) - a1(uu,0) + eps) ;
			MM=bb;
					}		
		   if (quanttype==2){	
			ecdf_temp[bb] = a1(uu,2) ;
			MM=bb;
					}				
			bb ++ ;
					}  else {
			uu ++ ;
					}
			} // end while
	// Rcpp::Rcout << "MM= " <<  MM <<  std::flush << std::endl ;		
	   if ( MM < BB  ){
		  hh=MM+1 ;
		  while (hh < BB ){   	  	  
		   ecdf_temp[bb] = a1(ngg-1,3) ;
		   hh ++ ;
			}
		   }
	
		for (int bb=0;bb<BB;bb++){
				ecdfMtemp[ bb + gg*BB + vv*GG*BB ] = ecdf_temp[bb] ;    	
			}  // end bb
			
	}  // end gg
	
	} // end vv
	
	return ( wrap( ecdfMtemp) ) ;
	}
	
	
//**************************************	
//**** logistic regression	
Rcpp::List bifie_estlogistic_helper( Rcpp::NumericVector y ,
	Rcpp::NumericMatrix X, Rcpp::NumericVector wgt ,
	Rcpp::NumericVector beta0 , double eps , int maxiter ){

int N=X.nrow();
int P=X.ncol();
double t1=0;

//*** create matrices in Armadillo
// design matrix X
arma::mat Xa0(N,P);
arma::mat Xa(N,P);
for (int nn=0;nn<N;nn++){
   for (int pp=0;pp<P;pp++){
	Xa0(nn,pp) = X(nn,pp) ;
			}
		 }
// outcome matrix y
arma::colvec ya(N) ;
for (int nn=0;nn<N;nn++){
	ya(nn,0) = y[nn] ;
		 }
// regression coefficients
arma::colvec beta_old(P);
arma::colvec beta_new(P);
for (int pp=0;pp<P;pp++){
	beta_old(pp,0) = beta0[pp] ;
			}

// temporary values in iterations
arma::colvec pred_logit(N);
arma::colvec ypred(N);
arma::colvec z(N);
arma::colvec AM(N);
arma::colvec wgta(N);
double pardiff=100;

int ii=0;

while( ( pardiff > eps ) & ( ii < maxiter ) ){

// within an iteration

// calculate predicted logit value and probability
for( int nn=0; nn <N; nn++){
pred_logit(nn,0)=0;
for ( int pp=0; pp <P; pp++){
	pred_logit(nn,0) += Xa0(nn,pp) * beta_old(pp,0) ;
		}
if ( pred_logit(nn,0) < - 15 ){
	pred_logit(nn,0) = - 15 ;
				}	
ypred(nn,0) = 1 / ( 1 + exp( - pred_logit(nn,0) ) ) ;
}
// calculate entries for A matrix and outcome z
for (int nn=0;nn<N;nn++){
	AM(nn,0) = ypred(nn,0) * ( 1 - ypred(nn,0) ) ;
	wgta(nn,0) = sqrt( AM(nn,0) * wgt[nn] ) ;
	z(nn,0) = pred_logit(nn,0) + ( ya(nn,0) - ypred(nn,0) )/AM(nn,0) ;
	z(nn,0) = wgta(nn,0) * z(nn,0) ;
		}
for (int nn=0;nn<N;nn++){
   for (int pp=0;pp<P;pp++){
   	   Xa(nn,pp)=Xa0(nn,pp)*wgta(nn,0);
   	   		}
   	   	}
// coefficient
beta_new = arma::solve(Xa, z);      // fit model y ~ X
// parameter difference
pardiff=0;
for (int pp=0;pp<P;pp++){
	t1 = beta_old(pp,0) - beta_new(pp,0) ;
	if (t1 < 0 ){ t1 = -t1 ; }
	if (t1 > pardiff){ pardiff = t1 ; }
		}
for (int pp=0;pp<P; pp++){
	beta_old(pp,0) = beta_new(pp,0);
		}
ii ++ ;		
	}
	
//*****************    
// OUTPUT            
        
return Rcpp::List::create( 
    _["pardiff"] = pardiff ,
    _["beta"] = beta_new ,
    _["iter"] = ii
    ) ;  
// maximal list length is 20!

	}
//*****************************************************

// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>

using namespace Rcpp;

// bifie_jack_timss
Rcpp::NumericMatrix bifie_jack_timss(Rcpp::NumericVector wgt, Rcpp::NumericVector jkzone, Rcpp::NumericVector jkrep, int RR, double jkfac, Rcpp::NumericVector prbar);
RcppExport SEXP BIFIEsurvey_bifie_jack_timss(SEXP wgtSEXP, SEXP jkzoneSEXP, SEXP jkrepSEXP, SEXP RRSEXP, SEXP jkfacSEXP, SEXP prbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wgt(wgtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type jkzone(jkzoneSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type jkrep(jkrepSEXP);
    Rcpp::traits::input_parameter< int >::type RR(RRSEXP);
    Rcpp::traits::input_parameter< double >::type jkfac(jkfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prbar(prbarSEXP);
    __result = Rcpp::wrap(bifie_jack_timss(wgt, jkzone, jkrep, RR, jkfac, prbar));
    return __result;
END_RCPP
}
// bifie_boot
Rcpp::List bifie_boot(Rcpp::NumericVector cumwgt, Rcpp::NumericMatrix rand_wgt);
RcppExport SEXP BIFIEsurvey_bifie_boot(SEXP cumwgtSEXP, SEXP rand_wgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cumwgt(cumwgtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rand_wgt(rand_wgtSEXP);
    __result = Rcpp::wrap(bifie_boot(cumwgt, rand_wgt));
    return __result;
END_RCPP
}
// bifie_bifiedata2bifiecdata
Rcpp::List bifie_bifiedata2bifiecdata(Rcpp::NumericMatrix datalistM, int Nimp);
RcppExport SEXP BIFIEsurvey_bifie_bifiedata2bifiecdata(SEXP datalistMSEXP, SEXP NimpSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalistM(datalistMSEXP);
    Rcpp::traits::input_parameter< int >::type Nimp(NimpSEXP);
    __result = Rcpp::wrap(bifie_bifiedata2bifiecdata(datalistM, Nimp));
    return __result;
END_RCPP
}
// bifie_bifiecdata2bifiedata
Rcpp::List bifie_bifiecdata2bifiedata(Rcpp::NumericMatrix datalistM_ind, Rcpp::NumericMatrix datalistM_imputed, int Nimp, Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix datalistM_impindex);
RcppExport SEXP BIFIEsurvey_bifie_bifiecdata2bifiedata(SEXP datalistM_indSEXP, SEXP datalistM_imputedSEXP, SEXP NimpSEXP, SEXP dat1SEXP, SEXP datalistM_impindexSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalistM_ind(datalistM_indSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalistM_imputed(datalistM_imputedSEXP);
    Rcpp::traits::input_parameter< int >::type Nimp(NimpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat1(dat1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalistM_impindex(datalistM_impindexSEXP);
    __result = Rcpp::wrap(bifie_bifiecdata2bifiedata(datalistM_ind, datalistM_imputed, Nimp, dat1, datalistM_impindex));
    return __result;
END_RCPP
}
// bifie_bifiedata_stepwise
Rcpp::List bifie_bifiedata_stepwise(Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix dat_ind, int Nmiss);
RcppExport SEXP BIFIEsurvey_bifie_bifiedata_stepwise(SEXP dat1SEXP, SEXP dat_indSEXP, SEXP NmissSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat1(dat1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat_ind(dat_indSEXP);
    Rcpp::traits::input_parameter< int >::type Nmiss(NmissSEXP);
    __result = Rcpp::wrap(bifie_bifiedata_stepwise(dat1, dat_ind, Nmiss));
    return __result;
END_RCPP
}
// univar_multiple_V2group
Rcpp::List univar_multiple_V2group(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values);
RcppExport SEXP BIFIEsurvey_univar_multiple_V2group(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    __result = Rcpp::wrap(univar_multiple_V2group(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values));
    return __result;
END_RCPP
}
// bifie_freq
Rcpp::List bifie_freq(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericMatrix vars_values, Rcpp::NumericVector vars_values_numb);
RcppExport SEXP BIFIEsurvey_bifie_freq(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP vars_valuesSEXP, SEXP vars_values_numbSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type vars_values(vars_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_values_numb(vars_values_numbSEXP);
    __result = Rcpp::wrap(bifie_freq(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, vars_values, vars_values_numb));
    return __result;
END_RCPP
}
// bifie_correl
Rcpp::List bifie_correl(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values);
RcppExport SEXP BIFIEsurvey_bifie_correl(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    __result = Rcpp::wrap(bifie_correl(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values));
    return __result;
END_RCPP
}
// bifie_linreg
Rcpp::List bifie_linreg(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector dep_index, Rcpp::NumericVector pre_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values);
RcppExport SEXP BIFIEsurvey_bifie_linreg(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP dep_indexSEXP, SEXP pre_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type dep_index(dep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pre_index(pre_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    __result = Rcpp::wrap(bifie_linreg(datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values));
    return __result;
END_RCPP
}
// bifie_waldtest
Rcpp::List bifie_waldtest(Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols, Rcpp::NumericVector fayfac);
RcppExport SEXP BIFIEsurvey_bifie_waldtest(SEXP parsMSEXP, SEXP parsrepMSEXP, SEXP CdesSEXP, SEXP rdesSEXP, SEXP CcolsSEXP, SEXP fayfacSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsM(parsMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsrepM(parsrepMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Cdes(CdesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rdes(rdesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Ccols(CcolsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    __result = Rcpp::wrap(bifie_waldtest(parsM, parsrepM, Cdes, rdes, Ccols, fayfac));
    return __result;
END_RCPP
}
// bifie_comp_vcov_within
Rcpp::List bifie_comp_vcov_within(Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, Rcpp::NumericVector fayfac, int RR, int Nimp);
RcppExport SEXP BIFIEsurvey_bifie_comp_vcov_within(SEXP parsMSEXP, SEXP parsrepMSEXP, SEXP fayfacSEXP, SEXP RRSEXP, SEXP NimpSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsM(parsMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsrepM(parsrepMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< int >::type RR(RRSEXP);
    Rcpp::traits::input_parameter< int >::type Nimp(NimpSEXP);
    __result = Rcpp::wrap(bifie_comp_vcov_within(parsM, parsrepM, fayfac, RR, Nimp));
    return __result;
END_RCPP
}
// bifie_comp_vcov
Rcpp::List bifie_comp_vcov(Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols, Rcpp::NumericVector fayfac);
RcppExport SEXP BIFIEsurvey_bifie_comp_vcov(SEXP parsMSEXP, SEXP parsrepMSEXP, SEXP CdesSEXP, SEXP rdesSEXP, SEXP CcolsSEXP, SEXP fayfacSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsM(parsMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type parsrepM(parsrepMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Cdes(CdesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rdes(rdesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Ccols(CcolsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    __result = Rcpp::wrap(bifie_comp_vcov(parsM, parsrepM, Cdes, rdes, Ccols, fayfac));
    return __result;
END_RCPP
}
// bifie_test_univar
Rcpp::List bifie_test_univar(Rcpp::NumericMatrix mean1M, Rcpp::NumericMatrix sd1M, Rcpp::NumericMatrix sumweightM, int GG, Rcpp::NumericVector group_values, Rcpp::NumericMatrix mean1repM, Rcpp::NumericMatrix sd1repM, Rcpp::NumericMatrix sumweightrepM, Rcpp::NumericVector fayfac);
RcppExport SEXP BIFIEsurvey_bifie_test_univar(SEXP mean1MSEXP, SEXP sd1MSEXP, SEXP sumweightMSEXP, SEXP GGSEXP, SEXP group_valuesSEXP, SEXP mean1repMSEXP, SEXP sd1repMSEXP, SEXP sumweightrepMSEXP, SEXP fayfacSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mean1M(mean1MSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type sd1M(sd1MSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type sumweightM(sumweightMSEXP);
    Rcpp::traits::input_parameter< int >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mean1repM(mean1repMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type sd1repM(sd1repMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type sumweightrepM(sumweightrepMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    __result = Rcpp::wrap(bifie_test_univar(mean1M, sd1M, sumweightM, GG, group_values, mean1repM, sd1repM, sumweightrepM, fayfac));
    return __result;
END_RCPP
}
// bifie_crosstab
Rcpp::List bifie_crosstab(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_values1, Rcpp::NumericVector vars_index1, Rcpp::NumericVector vars_values2, Rcpp::NumericVector vars_index2, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values);
RcppExport SEXP BIFIEsurvey_bifie_crosstab(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_values1SEXP, SEXP vars_index1SEXP, SEXP vars_values2SEXP, SEXP vars_index2SEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_values1(vars_values1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index1(vars_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_values2(vars_values2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index2(vars_index2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    __result = Rcpp::wrap(bifie_crosstab(datalist, wgt1, wgtrep, vars_values1, vars_index1, vars_values2, vars_index2, fayfac, NI, group_index1, group_values));
    return __result;
END_RCPP
}
// bifie_by
Rcpp::List bifie_by(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::Function userfct);
RcppExport SEXP BIFIEsurvey_bifie_by(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP userfctSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type userfct(userfctSEXP);
    __result = Rcpp::wrap(bifie_by(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, userfct));
    return __result;
END_RCPP
}
// bifie_hist
Rcpp::List bifie_hist(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericVector breaks);
RcppExport SEXP BIFIEsurvey_bifie_hist(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type breaks(breaksSEXP);
    __result = Rcpp::wrap(bifie_hist(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks));
    return __result;
END_RCPP
}
// bifie_ecdf
Rcpp::List bifie_ecdf(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericVector breaks, int quanttype, int maxval);
RcppExport SEXP BIFIEsurvey_bifie_ecdf(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP breaksSEXP, SEXP quanttypeSEXP, SEXP maxvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< int >::type quanttype(quanttypeSEXP);
    Rcpp::traits::input_parameter< int >::type maxval(maxvalSEXP);
    __result = Rcpp::wrap(bifie_ecdf(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks, quanttype, maxval));
    return __result;
END_RCPP
}
// bifie_logistreg
Rcpp::List bifie_logistreg(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector dep_index, Rcpp::NumericVector pre_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, double eps, int maxiter);
RcppExport SEXP BIFIEsurvey_bifie_logistreg(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP dep_indexSEXP, SEXP pre_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP epsSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type dep_index(dep_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pre_index(pre_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    __result = Rcpp::wrap(bifie_logistreg(datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values, eps, maxiter));
    return __result;
END_RCPP
}
// bifie_fasttable
Rcpp::List bifie_fasttable(Rcpp::NumericMatrix datavec);
RcppExport SEXP BIFIEsurvey_bifie_fasttable(SEXP datavecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datavec(datavecSEXP);
    __result = Rcpp::wrap(bifie_fasttable(datavec));
    return __result;
END_RCPP
}
// bifie_table1_character
Rcpp::List bifie_table1_character(Rcpp::CharacterVector datavec);
RcppExport SEXP BIFIEsurvey_bifie_table1_character(SEXP datavecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type datavec(datavecSEXP);
    __result = Rcpp::wrap(bifie_table1_character(datavec));
    return __result;
END_RCPP
}
// bifie_mla2
Rcpp::List bifie_mla2(Rcpp::NumericMatrix X_list, Rcpp::NumericMatrix Z_list, Rcpp::NumericVector y_list, Rcpp::NumericVector wgttot, Rcpp::NumericVector wgtlev2, Rcpp::NumericVector wgtlev1, double globconv, int maxiter, Rcpp::NumericVector group, Rcpp::NumericVector group_values, Rcpp::NumericVector cluster, Rcpp::NumericMatrix wgtrep, int Nimp, Rcpp::NumericVector fayfac, Rcpp::NumericMatrix recov_constraint, int is_rcov_constraint);
RcppExport SEXP BIFIEsurvey_bifie_mla2(SEXP X_listSEXP, SEXP Z_listSEXP, SEXP y_listSEXP, SEXP wgttotSEXP, SEXP wgtlev2SEXP, SEXP wgtlev1SEXP, SEXP globconvSEXP, SEXP maxiterSEXP, SEXP groupSEXP, SEXP group_valuesSEXP, SEXP clusterSEXP, SEXP wgtrepSEXP, SEXP NimpSEXP, SEXP fayfacSEXP, SEXP recov_constraintSEXP, SEXP is_rcov_constraintSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Z_list(Z_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_list(y_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wgttot(wgttotSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wgtlev2(wgtlev2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wgtlev1(wgtlev1SEXP);
    Rcpp::traits::input_parameter< double >::type globconv(globconvSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cluster(clusterSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< int >::type Nimp(NimpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type recov_constraint(recov_constraintSEXP);
    Rcpp::traits::input_parameter< int >::type is_rcov_constraint(is_rcov_constraintSEXP);
    __result = Rcpp::wrap(bifie_mla2(X_list, Z_list, y_list, wgttot, wgtlev2, wgtlev1, globconv, maxiter, group, group_values, cluster, wgtrep, Nimp, fayfac, recov_constraint, is_rcov_constraint));
    return __result;
END_RCPP
}
// bifie_pathmodel
Rcpp::List bifie_pathmodel(Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1, Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values, Rcpp::NumericMatrix L, Rcpp::NumericVector L_row_index, int NL, Rcpp::NumericMatrix E, Rcpp::NumericMatrix R, Rcpp::NumericVector R_row_index, Rcpp::NumericMatrix coeff_index, int NP0, Rcpp::NumericVector unreliability);
RcppExport SEXP BIFIEsurvey_bifie_pathmodel(SEXP datalistSEXP, SEXP wgt1SEXP, SEXP wgtrepSEXP, SEXP vars_indexSEXP, SEXP fayfacSEXP, SEXP NISEXP, SEXP group_index1SEXP, SEXP group_valuesSEXP, SEXP LSEXP, SEXP L_row_indexSEXP, SEXP NLSEXP, SEXP ESEXP, SEXP RSEXP, SEXP R_row_indexSEXP, SEXP coeff_indexSEXP, SEXP NP0SEXP, SEXP unreliabilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type datalist(datalistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgt1(wgt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wgtrep(wgtrepSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vars_index(vars_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fayfac(fayfacSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type NI(NISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_index1(group_index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type group_values(group_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type L_row_index(L_row_indexSEXP);
    Rcpp::traits::input_parameter< int >::type NL(NLSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type E(ESEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type R(RSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type R_row_index(R_row_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type coeff_index(coeff_indexSEXP);
    Rcpp::traits::input_parameter< int >::type NP0(NP0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type unreliability(unreliabilitySEXP);
    __result = Rcpp::wrap(bifie_pathmodel(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, L, L_row_index, NL, E, R, R_row_index, coeff_index, NP0, unreliability));
    return __result;
END_RCPP
}
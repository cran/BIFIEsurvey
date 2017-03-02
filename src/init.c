#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BIFIEsurvey_bifie_bifiecdata2bifiedata(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_bifiedata_stepwise(SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_bifiedata2bifiecdata(SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_boot(SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_by(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_comp_vcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_comp_vcov_within(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_correl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_crosstab(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_ecdf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_fasttable(SEXP);
extern SEXP BIFIEsurvey_bifie_freq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_hist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_jack_timss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_linreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_logistreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_mla2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_pathmodel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_table1_character(SEXP);
extern SEXP BIFIEsurvey_bifie_test_univar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_bifie_waldtest(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BIFIEsurvey_univar_multiple_V2group(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BIFIEsurvey_bifie_bifiecdata2bifiedata", (DL_FUNC) &BIFIEsurvey_bifie_bifiecdata2bifiedata,  5},
    {"BIFIEsurvey_bifie_bifiedata_stepwise",   (DL_FUNC) &BIFIEsurvey_bifie_bifiedata_stepwise,    3},
    {"BIFIEsurvey_bifie_bifiedata2bifiecdata", (DL_FUNC) &BIFIEsurvey_bifie_bifiedata2bifiecdata,  2},
    {"BIFIEsurvey_bifie_boot",                 (DL_FUNC) &BIFIEsurvey_bifie_boot,                  2},
    {"BIFIEsurvey_bifie_by",                   (DL_FUNC) &BIFIEsurvey_bifie_by,                    9},
    {"BIFIEsurvey_bifie_comp_vcov",            (DL_FUNC) &BIFIEsurvey_bifie_comp_vcov,             6},
    {"BIFIEsurvey_bifie_comp_vcov_within",     (DL_FUNC) &BIFIEsurvey_bifie_comp_vcov_within,      5},
    {"BIFIEsurvey_bifie_correl",               (DL_FUNC) &BIFIEsurvey_bifie_correl,                8},
    {"BIFIEsurvey_bifie_crosstab",             (DL_FUNC) &BIFIEsurvey_bifie_crosstab,             11},
    {"BIFIEsurvey_bifie_ecdf",                 (DL_FUNC) &BIFIEsurvey_bifie_ecdf,                 11},
    {"BIFIEsurvey_bifie_fasttable",            (DL_FUNC) &BIFIEsurvey_bifie_fasttable,             1},
    {"BIFIEsurvey_bifie_freq",                 (DL_FUNC) &BIFIEsurvey_bifie_freq,                 10},
    {"BIFIEsurvey_bifie_hist",                 (DL_FUNC) &BIFIEsurvey_bifie_hist,                  9},
    {"BIFIEsurvey_bifie_jack_timss",           (DL_FUNC) &BIFIEsurvey_bifie_jack_timss,            6},
    {"BIFIEsurvey_bifie_linreg",               (DL_FUNC) &BIFIEsurvey_bifie_linreg,                9},
    {"BIFIEsurvey_bifie_logistreg",            (DL_FUNC) &BIFIEsurvey_bifie_logistreg,            11},
    {"BIFIEsurvey_bifie_mla2",                 (DL_FUNC) &BIFIEsurvey_bifie_mla2,                 16},
    {"BIFIEsurvey_bifie_pathmodel",            (DL_FUNC) &BIFIEsurvey_bifie_pathmodel,            17},
    {"BIFIEsurvey_bifie_table1_character",     (DL_FUNC) &BIFIEsurvey_bifie_table1_character,      1},
    {"BIFIEsurvey_bifie_test_univar",          (DL_FUNC) &BIFIEsurvey_bifie_test_univar,           9},
    {"BIFIEsurvey_bifie_waldtest",             (DL_FUNC) &BIFIEsurvey_bifie_waldtest,              6},
    {"BIFIEsurvey_univar_multiple_V2group",    (DL_FUNC) &BIFIEsurvey_univar_multiple_V2group,     8},
    {NULL, NULL, 0}
};

void R_init_BIFIEsurvey(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

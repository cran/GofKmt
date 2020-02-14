// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// KmtMain
List KmtMain(arma::vec X, arma::mat NormalMat, arma::mat LogisMat, arma::mat ReMat, arma::mat CauchyMat, String strDist, int bGraph, int nNum);
RcppExport SEXP _GofKmt_KmtMain(SEXP XSEXP, SEXP NormalMatSEXP, SEXP LogisMatSEXP, SEXP ReMatSEXP, SEXP CauchyMatSEXP, SEXP strDistSEXP, SEXP bGraphSEXP, SEXP nNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type NormalMat(NormalMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LogisMat(LogisMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ReMat(ReMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type CauchyMat(CauchyMatSEXP);
    Rcpp::traits::input_parameter< String >::type strDist(strDistSEXP);
    Rcpp::traits::input_parameter< int >::type bGraph(bGraphSEXP);
    Rcpp::traits::input_parameter< int >::type nNum(nNumSEXP);
    rcpp_result_gen = Rcpp::wrap(KmtMain(X, NormalMat, LogisMat, ReMat, CauchyMat, strDist, bGraph, nNum));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GofKmt_KmtMain", (DL_FUNC) &_GofKmt_KmtMain, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_GofKmt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
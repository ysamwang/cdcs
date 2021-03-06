// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bnbHelperMulti
Rcpp::List bnbHelperMulti(const arma::uvec ancest, const arma::uvec children, const arma::mat& dat, arma::ivec K, int aggType, int bs, int intercept);
RcppExport SEXP _cdcs_bnbHelperMulti(SEXP ancestSEXP, SEXP childrenSEXP, SEXP datSEXP, SEXP KSEXP, SEXP aggTypeSEXP, SEXP bsSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec >::type ancest(ancestSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type children(childrenSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type aggType(aggTypeSEXP);
    Rcpp::traits::input_parameter< int >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(bnbHelperMulti(ancest, children, dat, K, aggType, bs, intercept));
    return rcpp_result_gen;
END_RCPP
}
// singleTestcppMulti
Rcpp::List singleTestcppMulti(const arma::mat& dat, const arma::colvec& Y, arma::ivec K, int aggType, int bs, int intercept);
RcppExport SEXP _cdcs_singleTestcppMulti(SEXP datSEXP, SEXP YSEXP, SEXP KSEXP, SEXP aggTypeSEXP, SEXP bsSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type aggType(aggTypeSEXP);
    Rcpp::traits::input_parameter< int >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(singleTestcppMulti(dat, Y, K, aggType, bs, intercept));
    return rcpp_result_gen;
END_RCPP
}
// exhaustiveHelperMulti
Rcpp::List exhaustiveHelperMulti(const arma::uvec ordered, const arma::mat& dat, arma::ivec K, int aggType, int bs, int intercept);
RcppExport SEXP _cdcs_exhaustiveHelperMulti(SEXP orderedSEXP, SEXP datSEXP, SEXP KSEXP, SEXP aggTypeSEXP, SEXP bsSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec >::type ordered(orderedSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type aggType(aggTypeSEXP);
    Rcpp::traits::input_parameter< int >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< int >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(exhaustiveHelperMulti(ordered, dat, K, aggType, bs, intercept));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cdcs_bnbHelperMulti", (DL_FUNC) &_cdcs_bnbHelperMulti, 7},
    {"_cdcs_singleTestcppMulti", (DL_FUNC) &_cdcs_singleTestcppMulti, 6},
    {"_cdcs_exhaustiveHelperMulti", (DL_FUNC) &_cdcs_exhaustiveHelperMulti, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_cdcs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

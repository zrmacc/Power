// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fitOLS
SEXP fitOLS(const arma::colvec y, const arma::mat X);
RcppExport SEXP _Power_fitOLS(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fitOLS(y, X));
    return rcpp_result_gen;
END_RCPP
}
// matIP
SEXP matIP(const arma::mat A, const arma::mat B);
RcppExport SEXP _Power_matIP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matIP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matQF
SEXP matQF(const arma::mat X, const arma::mat A);
RcppExport SEXP _Power_matQF(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matQF(X, A));
    return rcpp_result_gen;
END_RCPP
}
// matInv
SEXP matInv(const arma::mat A);
RcppExport SEXP _Power_matInv(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matInv(A));
    return rcpp_result_gen;
END_RCPP
}
// SchurC
SEXP SchurC(const arma::mat Ibb, const arma::mat Iaa, const arma::mat Iba);
RcppExport SEXP _Power_SchurC(SEXP IbbSEXP, SEXP IaaSEXP, SEXP IbaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Ibb(IbbSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Iaa(IaaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Iba(IbaSEXP);
    rcpp_result_gen = Rcpp::wrap(SchurC(Ibb, Iaa, Iba));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Power_fitOLS", (DL_FUNC) &_Power_fitOLS, 2},
    {"_Power_matIP", (DL_FUNC) &_Power_matIP, 2},
    {"_Power_matQF", (DL_FUNC) &_Power_matQF, 2},
    {"_Power_matInv", (DL_FUNC) &_Power_matInv, 1},
    {"_Power_SchurC", (DL_FUNC) &_Power_SchurC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Power(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

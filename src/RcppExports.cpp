// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// countCrossPoly
double countCrossPoly(int dim, int r);
RcppExport SEXP rrSample_countCrossPoly(SEXP dimSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    __result = Rcpp::wrap(countCrossPoly(dim, r));
    return __result;
END_RCPP
}
// countIntPoints
Rcpp::String countIntPoints(arma::mat constMat, std::vector<std::string> rhs);
RcppExport SEXP rrSample_countIntPoints(SEXP constMatSEXP, SEXP rhsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type rhs(rhsSEXP);
    __result = Rcpp::wrap(countIntPoints(constMat, rhs));
    return __result;
END_RCPP
}
// estimateDiam
int estimateDiam(IntegerVector initial, IntegerMatrix moves);
RcppExport SEXP rrSample_estimateDiam(SEXP initialSEXP, SEXP movesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type moves(movesSEXP);
    __result = Rcpp::wrap(estimateDiam(initial, moves));
    return __result;
END_RCPP
}
// estimateMixing
Rcpp::String estimateMixing(arma::uvec u, arma::mat constMat, arma::mat moves, int diam, std::string nIntPoints, double tol);
RcppExport SEXP rrSample_estimateMixing(SEXP uSEXP, SEXP constMatSEXP, SEXP movesSEXP, SEXP diamSEXP, SEXP nIntPointsSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::uvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type moves(movesSEXP);
    Rcpp::traits::input_parameter< int >::type diam(diamSEXP);
    Rcpp::traits::input_parameter< std::string >::type nIntPoints(nIntPointsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(estimateMixing(u, constMat, moves, diam, nIntPoints, tol));
    return __result;
END_RCPP
}
// fiberWalk
List fiberWalk(arma::uvec initial, arma::mat constMat, arma::mat moves, unsigned int diam, double length, bool showOutput);
RcppExport SEXP rrSample_fiberWalk(SEXP initialSEXP, SEXP constMatSEXP, SEXP movesSEXP, SEXP diamSEXP, SEXP lengthSEXP, SEXP showOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::uvec >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type constMat(constMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type moves(movesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type diam(diamSEXP);
    Rcpp::traits::input_parameter< double >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< bool >::type showOutput(showOutputSEXP);
    __result = Rcpp::wrap(fiberWalk(initial, constMat, moves, diam, length, showOutput));
    return __result;
END_RCPP
}
// sampleCrossPoly
IntegerVector sampleCrossPoly(int dim, int r, bool showOutput);
RcppExport SEXP rrSample_sampleCrossPoly(SEXP dimSEXP, SEXP rSEXP, SEXP showOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< bool >::type showOutput(showOutputSEXP);
    __result = Rcpp::wrap(sampleCrossPoly(dim, r, showOutput));
    return __result;
END_RCPP
}

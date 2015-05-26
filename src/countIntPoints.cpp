#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(algstat)]]
#include <string>
#include <Rdefines.h>

using namespace Rcpp;

// [[Rcpp::export]]
double countIntPoints(arma::mat constMat,arma::uvec rhs){

  Function getOption("getOption");
  Function count("count");
  double nIntPoints=0;

  SEXP opt=getOption("lattePath");
  //Rf_length checks if there is an element at index 0
  if(!Rf_isNull(opt) && Rf_length(opt) > 0) {
     //std::string opt=CHAR(STRING_ELT(Ropt,0));
     std::cout << "Count integer points with algstat not implemented yet" << std::endl;
     return 0;
     //use algstat.count here!
  } else {
     std::cout << "LattE is not loaded" << std::endl;
     return 0;
  }

   return nIntPoints;
}

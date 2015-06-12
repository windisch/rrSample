#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(algstat)]]
#include <string>
#include <vector>
#include <Rdefines.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::String countIntPoints(arma::mat constMat,std::vector<std::string> rhs){

  Function getOption("getOption");
  Function countFiber("countFiber");

  SEXP opt=getOption("lattePath");
  //Rf_length checks if there is an element at index 0
  if(!Rf_isNull(opt) && Rf_length(opt) > 0) {
    return as<std::string>(countFiber(constMat,rhs));
  } else {
     std::cout << "LattE is not loaded" << std::endl;
     return 0;
  }
}

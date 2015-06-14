#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(algstat)]]
#include <string>
#include <vector>
#include <Rdefines.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Rcpp::String countIntPoints(arma::mat constMat,SEXP rhs){

  Function getOption("getOption");
  Function countFiber("countFiber");

  SEXP opt=getOption("lattePath");
  //Rf_length checks if there is an element at index 0
  if(!Rf_isNull(opt) && Rf_length(opt) > 0) {
    SEXP numIntPoints=countFiber(constMat,rhs);
    switch( TYPEOF(numIntPoints) ) {
    case INTSXP: {
         return std::to_string(as<int>(numIntPoints));
    }
    case STRSXP: {
         return CHAR(STRING_ELT(numIntPoints,0)); 
         }
    default: {
         std::cout << "Unkown datatype from countFiber" << std::endl; 
         return 0;
         }
    }
   
  } else {
     std::cout << "LattE is not loaded" << std::endl;
     return 0;
  }
}

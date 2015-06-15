#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(algstat)]]
#include <string>
#include <vector>
#include <Rdefines.h>

using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//' @title Count integer points in a polytope
//'
//' @description This method counts the integer points in a polytope defined by the
//' matrix \eqn{A} and the right-hand side vector \eqn{b}. \code{countIntPoint} 
//' uses the method \code{countFiber} of the \code{algstat} package (which is
//' done there by using \code{LattE}). 
//' 
//' @param A matrix
//' @param b CharacterVector or NumericVector
//' @return the number of integer points in the polytope defined by
//' constMat and rhs. The number is returned in a character string.
//' @name countIntPoints
//' @seealso
//' \code{\link{countCrossPoly}}
//' \code{\link{countFiber}}
//' @usage
//' A<-matrix(c(1,1,1,1,1,0),2,3)
//' b<-c("14","20")
//' countIntPoints(A,b)
// [[Rcpp::export]]
Rcpp::String countIntPoints(arma::mat A,SEXP b){

  Function getOption("getOption");
  Function countFiber("countFiber");

  SEXP opt=getOption("lattePath");
  //Rf_length checks if there is an element at index 0
  if(!Rf_isNull(opt) && Rf_length(opt) > 0) {
    SEXP numIntPoints=countFiber(A,b);
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

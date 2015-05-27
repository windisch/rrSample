#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(algstat)]]
#include <string>
#include <Rdefines.h>
#include <gmp.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::String countIntPoints(arma::mat constMat,arma::uvec rhs){

    mpz_t nIntPoints;
    mpz_inits(nIntPoints,NULL);
      
  Function getOption("getOption");
  Function count("count");

  SEXP opt=getOption("lattePath");
  //Rf_length checks if there is an element at index 0
  if(!Rf_isNull(opt) && Rf_length(opt) > 0) {
     //std::string opt=CHAR(STRING_ELT(Ropt,0));
     std::cout << "Count integer points with algstat not implemented yet" << std::endl;
     return 0;

    mpz_set_str(nIntPoints, "342342341234", 10);
    //use here: mpz_set_str(nIntPoints,count(),10));
    std::string sIntPoints = mpz_get_str(NULL,10,nIntPoints);
    return sIntPoints;
  } else {
     std::cout << "LattE is not loaded" << std::endl;
     return 0;
  }

}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::uvec computeCellBounds(arma::uvec u,arma::mat constMat,std::string type="upper"){

arma::uvec cB;
if(type=="upper") {
   cB <<3600<<7200<<3600<<7200<<14400<<7200<<3600<<7200<<3600 << arma::endr;
   }

if(type=="lower") {
   cB <<0<<0<<0<<0<<0<<0<<0<<0<<0 << arma::endr;
   }

return cB;
}

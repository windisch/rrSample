#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::mat latticeComplement(arma::mat constMat){

arma::mat lC;
lC <<0<<0<<0<<0<<0<<0<<0<<0<<1 << arma::endr
<<0<<0<<0<<0<<0<<0<<0<<1<<-1 << arma::endr
<<0<<0<<0<<0<<0<<0<<1<<-1<<0 << arma::endr
<<0<<0<<0<<0<<0<<0<<-1<<0<<0 << arma::endr
<<0<<0<<0<<0<<0<<1<<0<<0<<-1 << arma::endr
<<0<<0<<0<<0<<1<<-1<<0<<-1<<1 << arma::endr
<<0<<0<<0<<1<<-1<<0<<-1<<1<<0 << arma::endr
<<0<<0<<0<<-1<<0<<0<<1<<0<<0 << arma::endr
<<0<<0<<1<<0<<0<<-1<<0<<0<<0 << arma::endr
<<0<<1<<-1<<0<<-1<<1<<0<<0<<0 << arma::endr
<<1<<-1<<0<<-1<<1<<0<<0<<0<<0 << arma::endr
<<-1<<0<<0<<1<<0<<0<<0<<0<<0 << arma::endr
<<0<<0<<-1<<0<<0<<0<<0<<0<<0 << arma::endr
<<0<<-1<<1<<0<<0<<0<<0<<0<<0 << arma::endr
<<-1<<1<<0<<0<<0<<0<<0<<0<<0 << arma::endr
<<1<<0<<0<<0<<0<<0<<0<<0<<0 << arma::endr;

return lC;
}

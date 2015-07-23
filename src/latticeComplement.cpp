#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::mat latticeComplement(arma::mat moves){
// This returns latticeComplements for some test cases



Function identical("identical");
arma::mat lC;
arma::mat M33;
arma::mat M44;

//indep-03-03 - basic moves
M33<<0<<0<<0<<0<<0<<1<<1<<1<<1<<arma::endr
<<0<<0<<0<<1<<1<<-1<<-1<<0<<0<<arma::endr
<<0<<0<<0<<-1<<-1<<0<<0<<-1<<-1<<arma::endr
<<0<<1<<1<<0<<0<<-1<<0<<-1<<0<<arma::endr
<<1<<-1<<0<<-1<<0<<1<<0<<0<<0<<arma::endr
<<-1<<0<<-1<<1<<0<<0<<0<<1<<0<<arma::endr 
<<0<<-1<<-1<<0<<0<<0<<-1<<0<<-1<<arma::endr 
<<-1<<1<<0<<0<<-1<<0<<1<<0<<0<<arma::endr 
<<1<<0<<1<<0<<1<<0<<0<<0<<1<<arma::endr;


//indep-04-04 - basic moves
M44<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<1<<1<<1<<1<<1<<1<<1<<1<<1<<arma::endr
<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<1<<1<<1<<1<<1<<1<<-1<<-1<<-1<<0<<0<<0<<0<<0<<0<<arma::endr
<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<1<<1<<1<<-1<<-1<<-1<<0<<0<<0<<0<<0<<0<<-1<<-1<<-1<<0<<0<<0<<arma::endr
<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<-1<<-1<<0<<0<<0<<-1<<-1<<-1<<0<<0<<0<<0<<0<<0<<-1<<-1<<-1<<arma::endr
<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<1<<1<<1<<1<<1<<1<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<-1<<0<<0<<arma::endr
<<0<<0<<0<<0<<0<<0<<0<<0<<1<<1<<1<<1<<-1<<-1<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<arma::endr
<<0<<0<<0<<0<<0<<0<<1<<1<<-1<<-1<<0<<0<<0<<0<<-1<<-1<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<0<<arma::endr
<<0<<0<<0<<0<<0<<0<<-1<<-1<<0<<0<<-1<<-1<<0<<0<<0<<0<<-1<<-1<<1<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<0<<0<<arma::endr
<<0<<0<<0<<1<<1<<1<<0<<0<<0<<0<<0<<0<<-1<<0<<-1<<0<<-1<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<-1<<0<<arma::endr
<<0<<1<<1<<-1<<0<<0<<0<<0<<-1<<0<<-1<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<arma::endr
<<1<<-1<<0<<0<<-1<<0<<-1<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<arma::endr
<<-1<<0<<-1<<0<<0<<-1<<1<<0<<0<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<0<<arma::endr
<<0<<0<<0<<-1<<-1<<-1<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<-1<<0<<-1<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<-1<<arma::endr
<<0<<-1<<-1<<1<<0<<0<<0<<0<<0<<-1<<0<<-1<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<0<<-1<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<arma::endr
<<-1<<1<<0<<0<<1<<0<<0<<-1<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<-1<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<0<<0<<0<<arma::endr
<<1<<0<<1<<0<<0<<1<<0<<1<<0<<0<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<1<<0<<0<<0<<0<<0<<1<<0<<0<<0<<0<<0<<0<<0<<0<<1<<arma::endr;

if(as<bool>(identical(moves,M44)))
{
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
}


if(as<bool>(identical(moves,M33)))
{

lC<<0<<0<<1<<1<<arma::endr
<<0<<1<<-1<<0<<arma::endr
<<0<<-1<<0<<-1<<arma::endr
<<1<<0<<-1<<0<<arma::endr
<<-1<<-1<<1<<0<<arma::endr
<<0<<1<<0<<0<<arma::endr
<<-1<<0<<0<<-1<<arma::endr
<<1<<0<<0<<0<<arma::endr
<<0<<0<<0<<1<<arma::endr;

}

return lC;
}
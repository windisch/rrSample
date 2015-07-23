#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::vec computeCellBounds(unsigned int diam,arma::mat adaptedMoves,std::string type="upper"){
// This returns upper and lower cell bounds for some test cases

Function identical("identical");

arma::mat AM33;
arma::mat AM44;
arma::vec cB;

//indep-03-03 - adapted basic moves
AM33<<0<<0<<1<<1<<arma::endr
<<0<<1<<-1<<0<<arma::endr
<<0<<-1<<0<<-1<<arma::endr
<<1<<0<<-1<<0<<arma::endr
<<-1<<-1<<1<<0<<arma::endr
<<0<<1<<0<<0<<arma::endr
<<-1<<0<<0<<-1<<arma::endr
<<1<<0<<0<<0<<arma::endr
<<0<<0<<0<<1<<arma::endr;

//indep-04-04 - adapted basic moves
AM44 <<0<<0<<0<<0<<0<<0<<0<<0<<1 << arma::endr
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

if(as<bool>(identical(adaptedMoves,AM44)))
{
   if(type=="upper") {
      cB <<36<<72<<36<<72<<144<<72<<36<<72<<36 << arma::endr;
      }
   if(type=="lower") {
      cB <<-36<<-72<<-36<<-72<<-144<<-72<<-36<<-72<<-36 << arma::endr;
      }
}


if(as<bool>(identical(adaptedMoves,AM33)))
{
   if(type=="upper") {
      cB <<9<<9<<18<<9<< arma::endr;
      }
   if(type=="lower") {
      cB <<-9<<-9<<-18<<-9<< arma::endr;
      }
}




//multiply cell bound with diameter
for(unsigned int i=0;i<cB.n_elem;i++) {
   cB[i]=diam*cB[i];    
    
}

return cB;
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double estimateMixing(arma::vec u,arma::mat moves,int diam){
//estimateMixing computes an upper bound on the mixing time
//diam should be an upper bound on diameter

  Function countCrossPoly("countCrossPoly");
  Function floor("floor");
  //double mixing;
  double nAdaptedMoves;
  
   if(arma::rank(moves)==moves.n_cols){
     //columns of moves are linear independent
     //in this case, the number of adapted moves coincides with the
     //number of elements in the corresponding cross poyltope
     nAdaptedMoves=as<double>(countCrossPoly(moves.n_cols,diam));
   }
   else
   {
     //columns of moves are linear dependent
     std::cout << "Linear dependent moves are not implemented yet" << std::endl; 
     return 0;
   }

  //nIntPoints=1; //Improve this by using count and --ehrhart-taylor!

  //mixing=nAdaptedMoves/(4*1);

  return nAdaptedMoves/4;
}

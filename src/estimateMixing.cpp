#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int estimateMixing(IntegerVector u,IntegerMatrix moves,int diam){
//estimateMixing computes an upper bound on the mxing time
//diam should be an upper bound on diameter
//TODO: Implement check on linear independence!

  Function countCrossPoly("countCrossPoly");
  Function floor("floor");
  int mixing=0;
  double nAdaptedMoves=0;
  double nIntPoints=0;
  
  nAdaptedMoves=as<double>(countCrossPoly(moves.ncol(),diam));
  nIntPoints=1;

  mixing=as<int>(floor(nAdaptedMoves/(4*nIntPoints)))+1;

  return mixing;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double estimateMixing(IntegerVector u,IntegerMatrix moves,int diam){
//estimateMixing computes an upper bound on the mxing time
//diam should be an upper bound on diameter
//TODO: Implement check on linear independence!

  Function countCrossPoly("countCrossPoly");
  Function floor("floor");
  //double mixing;
  double nAdaptedMoves;
  //double nIntPoints=1;
  
  nAdaptedMoves=as<double>(countCrossPoly(moves.ncol(),diam));

  //nIntPoints=1; //Improve this by using count and --ehrhart-taylor!

  //mixing=nAdaptedMoves/(4*1);

  return nAdaptedMoves/4;
}

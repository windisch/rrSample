#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <boost/math/special_functions/log1p.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
// [[Rcpp::depends(BH)]]
#include <gmpxx.h>
#include <string>

using namespace Rcpp;
namespace bmp = boost::multiprecision;
// [[Rcpp::export]]
Rcpp::String estimateMixing(arma::uvec u,arma::mat constMat,arma::mat moves,int diam,std::string nIntPoints="",double tol=0.25){
//estimateMixing computes an upper bound on the mixing time

  //check input
  if(u.n_elem!=moves.n_rows or u.n_elem!=constMat.n_cols){
     std::cout << "Wrong dimensions" << std::endl;
     return 0;
  }

  Function countCrossPoly("countCrossPoly");
  Function countIntPoints("countIntPoints");
  //double mixing;
  double nAdaptedMoves;
  arma::ivec rhs(constMat.n_rows);

  
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


  //compute right-hand side for computations in affine semigroup
   unsigned int k,j;
   #pragma omp parallel shared(constMat,rhs,u) private(k,j) 
   {
   #pragma omp for  schedule(static)
      for (k=0; k<rhs.n_elem; k++){
         rhs[k]=0.;
         for (j=0; j<constMat.n_cols; j++){
            rhs[k]=(rhs[k])+(constMat(k,j)*(u[j]));
         }
      }
   }

   //estimate integer points in (constMat,rhs))
   if(nIntPoints.size()==0){
       //convert rhs to List of strings!
       nIntPoints=as<std::string>(countIntPoints(constMat,rhs));
   }

  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > sIntPoints(nIntPoints);
  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > sAdaptedMoves(nAdaptedMoves);
  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > res;

  //boost::math::log1p(arg) computes log(arg+1)
  res= boost::math::log1p(tol/sqrt(sIntPoints)-1)/boost::math::log1p(-(sIntPoints*sIntPoints)/(8*sAdaptedMoves*sAdaptedMoves));
  return res.str();
}

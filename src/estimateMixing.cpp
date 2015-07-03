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
Rcpp::String estimateMixing(arma::uvec u,arma::mat constMat,arma::mat moves,int diam,std::string nAdaptedMoves="",std::string nIntPoints="",double tol=0.25,std::string type="TVD"){
//estimateMixing computes an upper bound on the mixing time

  //check input
  if(u.n_elem!=moves.n_rows or u.n_elem!=constMat.n_cols){
     std::cout << "Wrong dimensions" << std::endl;
     return 0;
  }

  Function countCrossPoly("countCrossPoly");
  if(arma::rank(moves)==moves.n_cols and nAdaptedMoves==""){
      double dAdap=as<double>(countCrossPoly(moves.n_cols,diam));
      nAdaptedMoves=std::to_string(dAdap);
  }

  //TODO: catch empty nAdaptedMoves for linear dependent moves

  Function countIntPoints("countIntPoints");

  mpz_class rhs[constMat.n_rows];

  //compute right-hand side for computations in affine semigroup
   unsigned int k,j;
   #pragma omp parallel shared(constMat,rhs,u) private(k,j) 
   {
   #pragma omp for  schedule(static)
      for (k=0; k<constMat.n_rows; k++){
         rhs[k]=mpz_class("0");
         for (j=0; j<constMat.n_cols; j++){
            rhs[k]=rhs[k]+mpz_class(constMat(k,j)*(u[j]));
         }
      }
   }

   //estimate integer points in (constMat,rhs))
   if(nIntPoints.size()==0){
       //convert rhs to Rcpp::List
       CharacterVector lrhs(constMat.n_rows);
       for (k=0; k<constMat.n_rows; k++){
           lrhs[k]=(rhs[k]).get_str();
       }

       nIntPoints=as<std::string>(countIntPoints(constMat,lrhs));
   }

  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > sIntPoints(nIntPoints);
  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > sAdaptedMoves(nAdaptedMoves);
  bmp::number<bmp::mpfr_float_backend<50,bmp::allocate_dynamic> > res;

  //boost::math::log1p(arg) computes log(arg+1)


  res= floor(boost::math::log1p(tol-1)/boost::math::log1p((sAdaptedMoves-sIntPoints)/sAdaptedMoves-1));
  return res.str();

  if(type=="TVD") {
  res= floor(boost::math::log1p(tol/sqrt(sIntPoints)-1)/boost::math::log1p(-(sIntPoints*sIntPoints)/(8*sAdaptedMoves*sAdaptedMoves)));
   }
   else {
       
  res= floor(boost::math::log1p(tol-1)/boost::math::log1p(-(sIntPoints*sIntPoints)/(8*sAdaptedMoves*sAdaptedMoves)));
       
       }
  return res.str();




}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <boost/math/common_factor.hpp>
// [[Rcpp::depends(BH)]]

using namespace Rcpp;

// [[Rcpp::export]]
double estimateMixing(arma::uvec u,arma::mat constMat,arma::mat moves,int diam){
//estimateMixing computes an upper bound on the mixing time
//diam should be an upper bound on diameter
//
   //check input
  if(u.n_elem!=moves.n_rows or u.n_elem!=constMat.n_cols){
     std::cout << "Wrong dimensions" << std::endl;
     return 0;
  }

  Function countCrossPoly("countCrossPoly");
  //double mixing;
  double nAdaptedMoves;
  int gcd;
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

   for (k=0; k<rhs.n_elem; k++){
      std::cout << rhs[k] << std::endl;
   }
   
   //get greatest common divisor of entries of rhs
   gcd=rhs[0];
   for (unsigned int i=1; i<rhs.n_elem; i++){
       gcd=boost::math::gcd(gcd,rhs[i]);
   }
      std::cout << "GCD:"  << gcd << std::endl;
   //right hand side to consider
   #pragma omp parallel for
   for (k=0; k<rhs.n_elem; k++){
      rhs[k]=rhs[k]/gcd;
   }

   std::cout << "NEW RHS" << std::endl;
   for (k=0; k<rhs.n_elem; k++){
      std::cout << rhs[k] << std::endl;
   }

   //try to estimate the integer points in (constMat,rhs)
   //this is a lower bound on int points in (constMat,gcd*rhs)



  //
  //
  //nIntPoints=1; //Improve this by using count and --ehrhart-taylor!

  //mixing=nAdaptedMoves/(4*1);

  return nAdaptedMoves/4;
}

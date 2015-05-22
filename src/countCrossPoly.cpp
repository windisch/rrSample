// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <boost/math/special_functions/binomial.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double countCrossPoly(int dim,int r){

   double count=0;
   double coeff=0;
   
   #pragma omp parallel for private(coeff)
   for(int k=0; k<=std::min(r,dim);k++){
       //TODO HERE: Catch overflows or write own functions within GMP
       coeff=boost::math::binomial_coefficient<double>(double(dim),double(k))*
       boost::math::binomial_coefficient<double>(double(r-k+dim),double(dim)); 
       #pragma omp critical 
       {
       count=count+coeff; 
       }
   }

  return count;
}

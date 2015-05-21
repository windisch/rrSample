// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <boost/math/special_functions/binomial.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sampleCrossPoly(int dim,int r){

  int remainder=r;
  DoubleVector probs(2*r+1);
  IntegerVector randomElement(dim);
  Function sample("sample");

  for(int i=0; i<dim, i++){

      //compute distribution on [-r+remainder,r-remainder]
      //
      //
      #pragma omp parallel for
      for(int j=0;j<probs.length();j++){
         probs[j]=0; 
      }

      #pragma omp parallel for
      for(int j=0;j<=remainder;j++)
      {
         countCrossPoly(d-1,r-remainder+j);     
         //write 1/() in probs at the right position!
      }

      //use sample function on [-remainder,remainder] and write output
      
      randomElement[i]=sample();
      remainder=remainder-randomElement[i];
      
  }

  return randomElement;
}

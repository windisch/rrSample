// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sampleCrossPoly(int dim,int r){

  int remainder=r;
  double prob=0;
  DoubleVector probs(2*r+1);
  IntegerVector states(2*r+1);
  IntegerVector randomElement(dim);
  Function sample("sample");
  Function countCrossPoly("countCrossPoly");

  for(int i=-r; i<=r; i++){
      states[i+r]=i;
  }

  for(int i=0; i<dim; i++){

      //initialize probability vector
      #pragma omp parallel for
      for(int j=0;j<probs.length();j++){
         probs[j]=0; 
      }

      #pragma omp parallel for private(prob)
      for(int j=0;j<=remainder;j++)
      {
         prob=as<double>(countCrossPoly(dim-i-1,r-remainder+j));     
         //write 1/() in probs at the right position!
         probs[r+i]=prob;
         probs[r-i]=prob;
      }
     
      for(int j=0;j<probs.length();j++){
         std::cout << i << "\t" << probs[j] << std::endl; 
      }


      //use sample function on [-remainder,remainder] and write output
      
      //randomElement[i]=sample();
      //remainder=remainder-randomElement[i];
      
  }

  return randomElement;
}

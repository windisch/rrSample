// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sampleCrossPoly(int dim,int r,bool showOutput=false){

  int remainder=r;
  double normConst=0;
  double prob=0;
  DoubleVector probs(2*r+1);
  IntegerVector states(2*r+1);
  IntegerVector randomElement(dim);
  Function sample("sample");
  Function countCrossPoly("countCrossPoly");

  #pragma omp parallel for
  for(int i=-r; i<=r; i++){
      states[i+r]=i;
  }

  for(int i=0; i<dim; i++){
   
      if(remainder>0){
         //initialize probability vector
         #pragma omp parallel for
         for(int j=0;j<probs.length();j++){
            probs[j]=0; 
         }
      
         normConst=0; 
         for(int j=0;j<=remainder;j++){
            prob=as<double>(countCrossPoly(dim-i-1,remainder-j));     
            probs[r+j]=prob;
            probs[r-j]=prob;
            if(j!=0){
               normConst=normConst+2*prob; 
            }
            else 
            {
               normConst=normConst+prob; 
            }

         }

         //compute normalizing constant
         #pragma omp parallel for
         for(int j=0;j<probs.length();j++){
            probs[j]=probs[j]/normConst;
         }

         if(showOutput){
            for(int j=0;j<probs.length();j++){
            std::cout << j << "\t" << probs[j] << std::endl; 
            }
            std::cout << "Normalizing Constant" << normConst << std::endl;
         }

         randomElement[i]=as<int>(sample(states,1,false,probs));

      } 
      else 
      {
         randomElement[i]=0; 
      }

      remainder=remainder-abs(randomElement[i]);
      if(showOutput){
         std::cout << "Random element" << "\t" << randomElement[i] << std::endl;
      }
  }

  return randomElement;
}

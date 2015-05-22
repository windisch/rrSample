#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

// [[Rcpp::export]]
List fiberWalk(arma::uvec initial,arma::mat constMat, arma::mat moves,unsigned int diam=0,double length=0,bool showOutput=false){


   if(arma::rank(moves)!=moves.n_cols){
     std::cout << "Linear independent moves needed" << std::endl; 
     return 0;
   }
   else
   {
     if(showOutput){
         std::cout << "Moves are linear independent" << std::endl; 
     }
   }

   //check input
   if(initial.n_elem!=moves.n_rows){
      std::cout << "Wrong dimensions" << std::endl;
      return 0;
   }

   //TODO: Implement  check on linear independence! Needed for now

  Function estimateDiam("estimateDiam");
  Function estimateMixing("estimateMixing");
  Function sampleCrossPoly("sampleCrossPoly");

  //estimate the diameter
  if(diam==0) {
      std::cout << "Estimate of diameter:";
      diam=as<unsigned int>(estimateDiam(initial,moves));
      std::cout << "\t" << diam << std::endl;
  }

  //estimate mixing
  if(length==0) {
     std::cout << "Estimate of mixing time:";
     length=as<double>(estimateMixing(initial,constMat,moves,diam));
     std::cout << "\t" << length << std::endl;
  }

  unsigned int dim = initial.n_elem;           // number of cells
  unsigned int N = moves.n_cols;               // number of moves
  int rejectionCounter=0;
  IntegerVector selection(1);
  IntegerVector proposal(dim);           
  IntegerVector current(dim);           
  bool applicable;
  IntegerVector move(dim);
  IntegerVector coeff(N);
 
  
  //this will fail if length flows integer values
  Progress p(length,true);

   #pragma omp parallel for
   for(unsigned int k=0;k<dim;k++){
      current[k]=initial[k];    
   }

  //this integer will overflow for large lengths
  for(double i = 0; i < length; ++i){
      
      if(!showOutput){
         p.increment();
      }

      //select move
      coeff=sampleCrossPoly(N,diam);

unsigned int k,j;
//matrix vector multiplication run in parallel
#pragma omp parallel shared(moves,move,coeff) private(k,j) 
{
#pragma omp for  schedule(static)
   for (k=0; k<dim; k++){
      move[k]=0.;
      for (j=0; j<N; j++){
         move[k]=(move[k])+(moves(k,j)*(coeff[j]));
      }
   }
}

      if(showOutput){
      std::cout << "Coefficient" << std::endl;
      for(unsigned int k=0; k<N; k++){
         std::cout << coeff[k] << "\t"; 
          }
          std::cout << std::endl;

      std::cout << "Move" << std::endl;
      for(unsigned int k=0; k<dim; k++){
         std::cout << move[k] << "\t";
          }
         std::cout << std::endl;
      }

      //check whether the move is applicable
      applicable = true;
      #pragma omp parallel for shared(applicable)
      for(unsigned int k = 0; k < dim; ++k){
        proposal[k] = current[k] + move[k];
        if(proposal[k]<0){
          #pragma omp critical
          {
           applicable=false;
          }
        }
      }

      //walk along edge or reject
      if(applicable){
          if(showOutput){
               std::cout << "traverse" << std::endl;
          }
      #pragma omp parallel for
        for(unsigned int k = 0; k < dim; ++k){
          current[k] = proposal[k];
        }
      } 
      else 
      {
          if(showOutput){
               std::cout << "reject" << std::endl;
          }
          rejectionCounter++;
      }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("sample") = current,
    Rcpp::Named("Rejections") = rejectionCounter,
    Rcpp::Named("Traversings") = length-rejectionCounter
  );

  return out;
}

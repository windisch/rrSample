// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
List fiberWalk(IntegerVector initial, IntegerMatrix moves,int diam, int length=0,bool showOutput=false){



  if(length==0) {
   // Implement an estimater of mixing time (volume of poyltope and
   // cross-poly
     length=1000; 
  }

  int dim = initial.size();             // number of cells
  int N = moves.ncol();               // number of moves
  int moveCounter=0;
  IntegerVector selection(1);
  IntegerVector proposal(dim);           
  IntegerVector current(dim);           
  //IntegerMatrix movesStack(dim,N);
  bool applicable;
  IntegerVector move(dim);
  IntegerVector coeff(N);
  Progress p(length,true);

  Function sampleCrossPoly("sampleCrossPoly");

   #pragma omp parallel for
   for(int k=0;k<dim;k++){
      current[k]=initial[k];    
   }


  for(int i = 0; i < length; ++i){

      p.increment();
      //std::cout << std::flush << i << "/" << length;

      //select move
      coeff=sampleCrossPoly(N,diam);

/*
     #pragma omp parallel for
     for(int j=0; j< N;j++) {
        for(int k=0; k<dim;k++) { 
        movesStack(k,j)=coeff[j]*moves(k,j);
        }
     }
     */

int k,j;
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
      for(int k=0; k<N; k++){
         std::cout << coeff[k] << "\t"; 
          }
          std::cout << std::endl;


      std::cout << "Move" << std::endl;
      for(int k=0; k<dim; k++){
         std::cout << move[k] << "\t";
          }
         std::cout << std::endl;
      }


      //check whether the move is applicable
      applicable = true;
      #pragma omp parallel for shared(applicable)
      for(int k = 0; k < dim; ++k){
        proposal[k] = current[k] + move[k];
        if(proposal[k]<0){
          #pragma omp critical
          {
           applicable=false;
          }
        }
      }

      //walk along edge
      if(applicable){
          moveCounter++;
      #pragma omp parallel for
        for(int k = 0; k < dim; ++k){
          current[k] = proposal[k];
        }
      }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("sample") = current,
    Rcpp::Named("NumMoves") = moveCounter,
    Rcpp::Named("Stucks") = length-moveCounter

    
  );

  return out;
}

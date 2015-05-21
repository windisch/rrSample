#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
List fiberWalk(IntegerVector initial, IntegerMatrix moves,int diam, int length=0){

  int dim = initial.size();             // number of cells
  int N = moves.ncol();               // number of moves
  IntegerVector selection(1);
  IntegerVector proposal(dim);           
  IntegerVector current(dim);           
  //IntegerMatrix movesStack(dim,N);
  bool applicable;
  IntegerVector move(dim);
  IntegerVector coeff(N);

  Function sampleCrossPoly("sampleCrossPoly");

   #pragma omp parallel for
   for(int i=0;i<dim;i++){
      current[i]=initial[i];    
   }

  if(length==0) {
   // Implement an estimater of mixing time (volume of poyltope and
   // cross-poly
     length=100; 
  }

  for(int i = 0; i < length; ++i){

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

      std::cout << "Coefficient" << std::endl;
      for(int k=0; k<N; k++){
         std::cout << coeff[k] << "\t"; 
          }
          std::cout << std::endl;


      std::cout << "Move" << std::endl;
      for(int k=0; k<dim; k++){
         std::cout << move[k] << std::endl;
          }


     /*std::cout << "Stack matrix" << std::endl;
     for(int k=0; k<dim;k++) { 
         for(int j=0; j< N;j++) {
            std::cout << movesStack(k,j) << "\t";
         }
         std::cout << std::endl;
       }
       */

      //compute proposal
      applicable = true;
      for(int k = 0; k < dim; ++k){
        proposal[k] = current[k] + move[k];
        if(proposal[k]<0){
           applicable=false;
           break;
        }
      }

      //walk along edge
      if(applicable){
      #pragma omp parallel for
        for(int k = 0; k < dim; ++k){
          current[k] = proposal[k];
        }
      }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("steps") = current
  );

  return out;
}

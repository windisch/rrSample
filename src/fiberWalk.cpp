#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]
#include <gmpxx.h>

using namespace Rcpp;

// [[Rcpp::export]]
List fiberWalk(arma::uvec initial,arma::mat constMat, arma::mat moves,unsigned int diam=0,SEXP length=R_NilValue,bool showOutput=false){


   //check moves on linear independence for now
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

  Function estimateDiam("estimateDiam");
  Function estimateMixing("estimateMixing");
  Function sampleCrossPoly("sampleCrossPoly");

  //estimate the diameter
  if(diam==0) {
      std::cout << "Estimation of diameter is not implemented yet" << std::endl;
      return 0;
      std::cout << "Estimate of diameter:";
      diam=as<unsigned int>(estimateDiam(initial,moves));
      std::cout << "\t" << diam << std::endl;
  }


  mpz_class steps;
  //estimate mixing
  if(length==R_NilValue) {
     SEXP estimate=estimateMixing(initial,constMat,moves,diam);
     steps=CHAR(STRING_ELT(estimate,0));

     std::cout << "Estimate of mixing time:";
     std::cout << "\t" << steps.get_str() << std::endl;
  } else {
      
    switch( TYPEOF(length) ) {
    case REALSXP: {
        steps=as<double>(length);
        break;
    }
    case INTSXP: {
         steps=std::to_string(as<int>(length));
         break;
    }
    case STRSXP: {
         steps=CHAR(STRING_ELT(length,0));
         break;
         }
    default: {
         std::cout << "Unkown datatype of length" << std::endl; 
         return 0;
         }
    }
  }

  unsigned int dim = initial.n_elem;           // number of cells
  unsigned int N = moves.n_cols;               // number of moves
  IntegerVector selection(1);
  IntegerVector proposal(dim);           
  IntegerVector current(dim);           
  bool applicable;
  IntegerVector move(dim);
  IntegerVector coeff(N);
 
  //this will fail if length flows integer values
  //Progress p(length,true);

   #pragma omp parallel for
   for(unsigned int k=0;k<dim;k++){
      current[k]=initial[k];    
   }

  mpz_class i("0",10);
  mpz_class rejectionCounter("0",10);
  mpz_class transitionCounter("0",10);

  while(i<steps){
      i=i+1;
      
   /*   if(!showOutput){
         p.increment();
      }*/

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
           //break the loop can speed things up in large dimension
          }
        }
      }

      //walk along edge or reject
      if(applicable){
          if(showOutput){
               std::cout << "traverse" << std::endl;
          }
          transitionCounter=transitionCounter+1;
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
          rejectionCounter=rejectionCounter+1;
      }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("sample") = current,
    Rcpp::Named("Rejections") = rejectionCounter.get_str(),
    //FIXME: insert steps here
    Rcpp::Named("Transitions") = transitionCounter.get_str()
  );

  return out;
}

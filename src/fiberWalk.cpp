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

   //check input
   if(initial.n_elem!=moves.n_rows){
      std::cout << "Wrong dimensions" << std::endl;
      return 0;
   }

  Function estimateDiam("estimateDiam");
  Function estimateMixing("estimateMixing");
  Function latticeComplement("latticeComplement");
  Function computeCellBounds("computeCellBounds");
  Function sampleCrossPoly("sampleCrossPoly");
  Function sample("sample");

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

  unsigned int dim = initial.n_elem;         
  IntegerVector proposal(dim);           
  IntegerVector current(dim);           
  bool applicable;
  IntegerVector move(dim);

  arma::mat adaptedMoves;
  arma::vec lower;
  arma::vec upper;
  unsigned int rank=0;
  bool linIndep;



  //check linear independence
  if(arma::rank(moves)!=moves.n_cols){
       linIndep=false;
   } else {
      linIndep=true;    
   }

   if(linIndep==false){
      //linear dependent moves; sample via lattice complement

      //compute Lattice complement
      //right now, this returs the PRODUCT M*L and not just L
      adaptedMoves = as<arma::mat>(latticeComplement(constMat));
      rank=adaptedMoves.n_cols;

      //compute cell bounds
      lower=as<arma::vec>(computeCellBounds(initial,constMat,"lower"));
      upper=as<arma::vec>(computeCellBounds(initial,constMat,"upper"));

   }
   else
   {
       //linear indepentent moves, sample from moves-matrix directly
       adaptedMoves=moves;
       rank=adaptedMoves.n_cols;
   }




  //initialize coefficients of appropriate size
  arma::vec coeff(rank);
  IntegerVector crossPolySample(rank);

  //this will fail if length flows integer values
  //Progress p(length,true);

   #pragma omp parallel for
   for(unsigned int k=0;k<dim;k++){
      current[k]=initial[k];    
   }

  mpz_class i("0",10);
  mpz_class rejectionCounter("0",10);
  mpz_class transitionCounter("0",10);

unsigned int k,j;

while(i<steps){
      i=i+1;

//sample coefficients
   if(linIndep==false) {
      //sample coefficients within cell bounds
      for(k=0;k<rank;k++) {
         coeff[k]=lower[k]+as<unsigned int>(sample(upper[k]-lower[k],1));
      }
   }
   else {
      //sample coefficients from cross-polytope
      crossPolySample=sampleCrossPoly(rank,diam);

      for(k=0;k<rank;k++) {
         coeff[k]=crossPolySample[k];
      }

   }

//compute move and apply it
#pragma omp parallel shared(adaptedMoves,move,coeff) private(k,j) 
{
#pragma omp for  schedule(static)
   for (k=0; k<dim; k++){
      move[k]=0.;
      for (j=0; j<rank; j++){
         move[k]=(move[k])+(adaptedMoves(k,j)*(coeff[j]));
      }
   }
}

      if(showOutput){
      std::cout << "Coefficient" << std::endl;
      for(k=0; k<rank; k++){
         std::cout << coeff[k] << "\t"; 
          }
          std::cout << std::endl;

      std::cout << "Move" << std::endl;
      for(k=0; k<dim; k++){
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

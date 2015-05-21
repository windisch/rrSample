#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
List fiberWalk(IntegerVector current, IntegerMatrix moves, int length){

  int n = current.size();             // number of cells
  int N = moves.ncol();               // number of moves
  IntegerVector selection(1);
  IntegerMatrix steps(n, length);         
  IntegerVector proposal(n);           
  bool applicable;
  IntegerVector move(n);

  Function sample("sample");


 // TODOS
 // Implement an estimater of mixing time (volume of poyltope and
 // volume of crossPoly)

  for(int i = 0; i < length; ++i){

      //select move
      std::cout << "select move" << std::endl;
      selection = sample(N, 1);
      std::cout << selection[0] << std::endl;

      //compute proposal
      applicable = true;
      for(int k = 0; k < n; ++k){
        proposal[k] = current[k] + moves(k, selection[0]-1);
        if(proposal[k]<0){
           applicable=false;
           break;
        }
      }

      //walk along edge
      if(applicable){
      #pragma omp parallel for
        for(int k = 0; k < n; ++k){
          current[k] = proposal[k];
        }
      }

    // assign state move
     #pragma omp parallel for
    for(int k = 0; k < n; ++k){
      steps(k,i) = current[k];
    }
  }

  // create out list
  List out = List::create(
    Rcpp::Named("steps") = current
  );

  return out;
}

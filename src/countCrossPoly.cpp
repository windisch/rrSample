// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <boost/math/special_functions/binomial.hpp>

using namespace Rcpp;

//' @title Count integer points in the cross polytope
//'
//' @description This method counts the integer points in 
//' the \eqn{r}th dilitation
//' of the \eqn{n}th dimensional cross-polytope. 
//' \deqn{ \{v\in Z^n: \sum_{i=1}^n|v_i|\le r \}}
//' This is done by evaluation the following formula 
//' \deqn{\sum_{k=0}^n\binom{n}{k}\cdot\binom{r-k+n}{n}}
//' which is exactly the cardinality of the set given above.
//'
//' @param dim integer
//' @param r integer 
//' @return the number of integer points in the cross polytope 
//' @name countCrossPoly
//' @references
//' \url{https://en.wikipedia.org/wiki/Cross-polytope}
//' @seealso
//' \code{\link{countIntPoints}}
//' \code{\link{countFiber}}
//' @usage
//' A<-matrix(c(1,1,1,1,1,0),2,3)
//' b<-c("14","20")
//' countIntPoints(A,b)
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

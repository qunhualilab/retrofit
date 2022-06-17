#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector retrofit_step2_rgamma(NumericVector shapes, NumericVector rates) {
/* Equivalent code
 */
  NumericVector ret (shapes.length());
  for(int i=0; i<shapes.length(); ++i){
    
    ret[i] = rgamma(1, shapes[i], 1/rates[i])[0];
  }
  
  return ret;
}

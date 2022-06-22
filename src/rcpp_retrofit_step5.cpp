#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector retrofit_step5_parameter_estimation(NumericVector original, 
                                                  NumericVector update, 
                                                  float rho) {
  
  NumericVector ret (original.length());
  for(int i=0; i<original.length(); ++i){
    ret[i] = (1-rho)*original[i] + rho*update[i];
  }
  return ret;
}

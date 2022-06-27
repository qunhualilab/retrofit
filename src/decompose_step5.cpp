#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void decompose_step5(NumericVector main, 
                     NumericVector asterisk, 
                     double rho) {
  
  for(int i=0; i<main.length(); ++i){
    main[i] = (1-rho)*main[i] + rho*asterisk[i];
  }
}

//sum.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_second_dim_sum(NumericVector v){
  NumericVector dim = v.attr("dim");
  int d1 = dim[0];
  int d2 = dim[1];
  int d3 = dim[2];
  
  int idx = 0;
  NumericVector sums (d1*d3);
  // double sum = 0;
  for(int k=0; k<d3; ++k){
    for(int j=0; j<d2; ++j){
      for(int i=0; i<d1; ++i){
        // Rprintf("[%i, %i, %i]: %f \n", i, j, k, v[idx]);
        sums[(k*d1+i)] += v[idx];
        // sum += v[idx];
        idx++;
      }
    }
  }
  return(sums);
}
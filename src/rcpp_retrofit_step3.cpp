#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector alpha_denominator(NumericVector phi_a_gks) {
/* Equivalent code 
 * for(s in 1:S){
 *  for(v in 1:G){
 *    phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
 *  }
 *}
 */
  //3d array
  NumericVector dim = phi_a_gks.attr("dim");
  int d1 = dim[0];
  int d2 = dim[1];
  int d3 = dim[2];
  
  //Calculate second dimension sums
  int idx = 0;
  NumericVector sums (d1*d3);
  for(int k=0; k<d3; ++k){
    for(int j=0; j<d2; ++j){
      for(int i=0; i<d1; ++i){
        sums[(k*d1+i)] += phi_a_gks[idx];
        idx++;
      }
    }
  }
  
  idx = 0;
  for(int k=0; k<d3; ++k){
    for(int j=0; j<d2; ++j){
      for(int i=0; i<d1; ++i){
        phi_a_gks[idx] = phi_a_gks[idx]/sums[(k*d1+i)];
        idx++;
      }
    }
  }
  
  return phi_a_gks;
}


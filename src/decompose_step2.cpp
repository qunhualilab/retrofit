#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
void decompose_step2(NumericVector shapes, 
                     NumericVector rates,
                     NumericVector &out) {
/* Equivalent code
 * for(k in 1:L){
 *    for(s in 1:S){
 *       H1[k,s]=rgamma(1, shape=eta_h[k,s], rate=gam_h[k,s])
 *    }
 *  }
 */
  for(int i=0; i<shapes.length(); ++i){
    out[i] = rgamma(1, shapes[i], 1/rates[i])[0];
  }
}

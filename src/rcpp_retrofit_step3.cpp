#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector phi_alpha_numerator(NumericVector W_gk, NumericVector TH_k, NumericVector H_ks, float lambda) {
  /* Equivalent code
   * for(s in 1:S){
   *  for(k in 1:K){
   *    phi_a_gks[,k,s]=((W_gk[,k] * TH_k[k]) +lamda)* H_ks[k,s]
   *  }
   * }
   */
  //dimensions
  // W_gk:(G,K), TH_K:(K), H_ks:(K,S)
  NumericVector dim_w = W_gk.attr("dim");
  int G = dim_w[0];
  int K = dim_w[1];
  NumericVector dim_h = H_ks.attr("dim");
  int S = dim_h[1];


  //Calculate numerator part
  int idx = 0;
  NumericVector phi_a_gks (G*K*S);
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        // w = W_gk[(k*G+g)];
        // th = TH_k[k];
        // h = H_ks[(s*K+k)];
        // phi_a_gks[idx] = (w*th + lambda)*h;
        phi_a_gks[idx] = (W_gk[(k*G+g)]*TH_k[k] + lambda)*H_ks[(s*K+k)];
        idx++;
      }
    }
  }

  return phi_a_gks;
}

// [[Rcpp::export]]
NumericVector phi_alpha_denominator(NumericVector phi_a_gks) {
/* Equivalent code
 * for(s in 1:S){
 *  for(v in 1:G){
 *    phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
 *  }
 *}
 */
  //3d array
  NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];

  //Calculate second dimension sums
  int idx = 0;
  NumericVector sums (G*S);
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        sums[(s*G+g)] += phi_a_gks[idx];
        idx++;
      }
    }
  }

  //Calculate denominator part
  idx = 0;
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        phi_a_gks[idx] = phi_a_gks[idx]/sums[(s*G+g)];
        idx++;
      }
    }
  }

  return phi_a_gks;
}


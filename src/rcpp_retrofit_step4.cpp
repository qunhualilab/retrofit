#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector retrofit_step4_alpha_w(NumericVector x_gs, NumericVector phi_a_gks, NumericVector phi_b_gk, float alpha_w_0) {
  /* Equivalent code
   * for(k in 1:K){
   *  alpha_W_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k]
   * }
   */
  //dimensions
  NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];

  //Calculate rowSums
  NumericVector alpha_w (G*K);
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        alpha_w[k*G+g] += x_gs[s*G+g]*phi_a_gks[s*G*K+k*G+g]*phi_b_gk[k*G+g];
      }
    }
  }
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      alpha_w[k*G+g] += alpha_w_0;
    }
  }

  return alpha_w;
}

// [[Rcpp::export]]
NumericVector retrofit_step4_alpha_h(NumericVector x_gs, NumericVector phi_a_gks, float alpha_h_0) {
  /* Equivalent code
   * for(k in 1:K){
   *  alpha_H_0 + colSums(x*phi_a_gks[,k,])
   * }
   */
  //dimensions
  NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];
  
  //Calculate colSums
  NumericVector alpha_h (K*S);
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        alpha_h[s*K+k] += x_gs[s*G+g]*phi_a_gks[s*G*K+k*G+g];
      }
    }
  }
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      alpha_h[s*K+k] += alpha_h_0;
    }
  }
  
  return alpha_h;
}

// [[Rcpp::export]]
NumericVector retrofit_step4_alpha_th(NumericVector x_gs, NumericVector phi_a_gks, NumericVector phi_b_gk, float alpha_th_0) {
  /* Equivalent code
   * for(k in 1:K){
   *  alpha_TH_0 + sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
   * }
   */
  //dimensions
  NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];
  
  //Calculate 
  NumericVector alpha_th (K);
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        alpha_th[k] += x_gs[s*G+g]*phi_a_gks[s*G*K+k*G+g]*phi_b_gk[k*G+g];
      }
    }
  }
  for(int k=0; k<K; ++k){
    alpha_th[k] += alpha_th_0;
  }
  
  return alpha_th;
}

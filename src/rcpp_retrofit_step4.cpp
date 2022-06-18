#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List retrofit_step4_alpha_calculation(NumericVector x_gs, 
                                      NumericVector phi_a_gks, 
                                      NumericVector phi_b_gk, 
                                      float alpha_w_0,
                                      float alpha_h_0,
                                      float alpha_th_0) {
  /* Equivalent code
   * for(k in 1:K){
   *  alpha_W_gk[,k]= alpha_W_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
   *  alpha_H_ks[k,]= alpha_H_0 + colSums(x*phi_a_gks[,k,])
   *  alpha_TH_k[k]= alpha_TH_0 + sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
   * }
   */
  //dimensions
  NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];

  //Calculate alpha_w
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
  
  //Calculate alpha_h
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
  
  //Calculate alpha_th
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
  

  List ret;
  ret["w"] = alpha_w;
  ret["h"] = alpha_h;
  ret["t"] = alpha_th;
  return ret;
}

// [[Rcpp::export]]
List retrofit_step4_beta_calculation(NumericVector W_gk,
                                     NumericVector H_ks,
                                     NumericVector TH_k,
                                     float beta_w_0,
                                     float beta_h_0,
                                     float beta_th_0,
                                     float lambda)
{
  return List();
}

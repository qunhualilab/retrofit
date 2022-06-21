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
  NumericVector alpha_h (K*S);
  NumericVector alpha_th (K);
  float x_gs_val;
  float phi_a_gks_val;
  float phi_b_gk_val;
  NumericVector::iterator phi_a_gks_iter = phi_a_gks.begin();
  NumericVector::iterator phi_b_gk_iter = phi_b_gk.begin();
  NumericVector::iterator x_gs_iter = x_gs.begin();
  NumericVector::iterator alpha_w_iter = alpha_w.begin();
  NumericVector::iterator alpha_h_iter = alpha_h.begin();
  NumericVector::iterator alpha_th_iter = alpha_th.begin();
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        // equivalent: x_gs_val = x_gs[s*G+g];
        x_gs_val = *x_gs_iter; ++x_gs_iter; if(g==G-1 && k!=K-1) {x_gs_iter -= G;}

        // equivalent: phi_a_gks_val = phi_a_gks[s*G*K+k*G+g];
        phi_a_gks_val = *phi_a_gks_iter; ++phi_a_gks_iter;

        // equivalent: phi_b_gk_val = phi_b_gk[k*G+g];
        phi_b_gk_val = *phi_b_gk_iter; ++phi_b_gk_iter; if(g==G-1 && k==K-1) {phi_b_gk_iter -= K*G;}

        // equivalent: alpha_w[k*G+g] += x_gs_val*phi_a_gks_val*phi_b_gk_val;
        *alpha_w_iter += x_gs_val*phi_a_gks_val*phi_b_gk_val;
         ++alpha_w_iter; if(g==G-1 && k==K-1) {alpha_w_iter -= K*G;}

        // equivalent: alpha_h[s*K+k] += x_gs_val*phi_a_gks_val;
        *alpha_h_iter += x_gs_val*phi_a_gks_val;
        if(g==G-1) {++alpha_h_iter;}

        // equivalent: alpha_th[k] += x_gs_val*phi_a_gks_val*phi_b_gk_val;
        *alpha_th_iter += x_gs_val*phi_a_gks_val*phi_b_gk_val;
        if(g==G-1) {++alpha_th_iter;} if(g==G-1 && k==K-1) {alpha_th_iter -= K;}
      }
    }
  }

  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      alpha_w[k*G+g] += alpha_w_0;
    }
  }
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      alpha_h[s*K+k] += alpha_h_0;
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
  /* Equivalent code
   * for(k in 1:K){
   * beta_W_gk[,k]= beta_W_0 + sum(H_ks[k,]*TH_k[k])
   * beta_H_ks[k,]= beta_H_0 + sum(W_gk[,k]*TH_k[k] + lamda)
   * beta_TH_k[k]= beta_TH_0 + sum(as.matrix(W_gk[,k]) %*% t(as.matrix(H_ks[k,])))
   * }
   */
  //dimensions
  NumericVector dim_w = W_gk.attr("dim");
  int G = dim_w[0];
  int K = dim_w[1];
  NumericVector dim_h = H_ks.attr("dim");
  int S = dim_h[1];
  
  NumericVector beta_w (G*K);
  NumericVector beta_h (K*S);
  NumericVector beta_th (K);
  NumericVector::iterator beta_w_iter = beta_w.begin();
  NumericVector::iterator beta_h_iter = beta_h.begin();
  NumericVector::iterator beta_th_iter = beta_th.begin();
  NumericVector::iterator W_gk_iter = W_gk.begin();
  NumericVector::iterator H_ks_iter = H_ks.begin();
  NumericVector::iterator TH_k_iter = TH_k.begin();
  double sum;
  
  // calculate beta_w
  // equivalent: beta_W_gk[,k]= beta_W_0 + sum(H_ks[k,]*TH_k[k])
  H_ks_iter = H_ks.begin();
  TH_k_iter = TH_k.begin();
  for(int k=0; k<K; ++k){

    sum = 0;
    for(int s=0; s<S; ++s){
      // equivalent: beta_w[k*G+g] += TH_k_val*H_ks_val;
      sum += (*TH_k_iter)*(*H_ks_iter);
      // equivalent: H_ks_val = H_ks[s*K+k];
      H_ks_iter += K; if(s==S-1){H_ks_iter -= K*S;}
    }
    ++H_ks_iter;
    ++TH_k_iter;

    for(int g=0; g<G; ++g){
      *beta_w_iter = sum;
      ++beta_w_iter;
    }
  }

  // calculate beta_h
  // equivalent: beta_H_ks[k,]= beta_H_0 + sum(W_gk[,k]*TH_k[k] + lamda)
  W_gk_iter = W_gk.begin();
  TH_k_iter = TH_k.begin();
  for(int k=0; k<K; ++k){
    sum = 0;

    for(int g=0; g<G; ++g){
      sum += (*W_gk_iter)*(*TH_k_iter)+lambda;
      // equivalent: W_gk_val = W_gk[k*G+g];
      ++W_gk_iter;
    }
    // equivalent: TH_k_val = TH_k[k];
    ++TH_k_iter;

    for(int s=0; s<S; ++s){
      *beta_h_iter = sum;
      beta_h_iter += K; if(s==S-1){beta_h_iter-=K*S;}
    }
    ++beta_h_iter;
  }


  W_gk_iter = W_gk.begin();
  H_ks_iter = H_ks.begin();
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        // equivalent: beta_th[k] += W_gk_val*H_ks_val;
        *beta_th_iter += (*W_gk_iter)*(*H_ks_iter);
        ++W_gk_iter;
      }
      ++H_ks_iter;
      ++beta_th_iter;
      if(k==K-1){beta_th_iter-=K;}
    }
    W_gk_iter -= G*K;
  }
    
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      beta_w[k*G+g] += beta_w_0;
    }
  }
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      beta_h[s*K+k] += beta_h_0;
    }
  }
  for(int k=0; k<K; ++k){
    beta_th[k] += beta_th_0;
  }
  
  
  List ret;
  ret["w"] = beta_w;
  ret["h"] = beta_h;
  ret["t"] = beta_th;
  return ret;
}

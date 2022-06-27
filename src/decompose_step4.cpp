#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List decompose_step4_alpha(NumericVector x_gs, 
                           List probabilities, 
                           double alpha_w_0,
                           double alpha_h_0,
                           double alpha_th_0,
                           NumericVector dim) {
  /* Equivalent code
   * for(k in 1:K){
   *  alpha_w_gk[,k]= alpha_W_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
   *  alpha_h_ks[k,]= alpha_H_0 + colSums(x*phi_a_gks[,k,])
   *  alpha_th_k[k]= alpha_TH_0 + sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
   * }
   */
  //dimensions
  // NumericVector dim = phi_a_gks.attr("dim");
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];

  //Calculate alpha_w
  NumericVector alpha_w (G*K);
  NumericVector alpha_h (K*S);
  NumericVector alpha_th (K);
  double x_gs_val;
  double phi_a_gks_val;
  double phi_b_gk_val;
  NumericVector phi_a_gks = probabilities["phi_a_gks"];
  NumericVector phi_b_gk = probabilities["phi_b_gk"];
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
List decompose_step4_beta(List distributions,
                          double beta_w_0,
                          double beta_h_0,
                          double beta_th_0,
                          double lambda,
                          NumericVector dim)
{
  /* Equivalent code
   * for(k in 1:K){
   * beta_w_gk[,k]= beta_W_0 + sum(h_ks[k,]*th_k[k])
   * beta_h_ks[k,]= beta_H_0 + sum(w_gk[,k]*th_k[k] + lamda)
   * beta_th_k[k]= beta_TH_0 + sum(as.matrix(w_gk[,k]) %*% t(as.matrix(h_ks[k,])))
   * }
   */
  //dimensions
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];
  
  NumericVector beta_w (G*K);
  NumericVector beta_h (K*S);
  NumericVector beta_th (K);
  NumericVector::iterator beta_w_iter = beta_w.begin();
  NumericVector::iterator beta_h_iter = beta_h.begin();
  NumericVector::iterator beta_th_iter = beta_th.begin();
  NumericVector w_gk = distributions["w_gk"];
  NumericVector h_ks = distributions["h_ks"];
  NumericVector th_k = distributions["th_k"];
  NumericVector::iterator w_gk_iter = w_gk.begin();
  NumericVector::iterator h_ks_iter = h_ks.begin();
  NumericVector::iterator th_k_iter = th_k.begin();
  double sum;
  
  // calculate beta_w
  // equivalent: beta_w_gk[,k]= beta_W_0 + sum(h_ks[k,]*th_k[k])
  h_ks_iter = h_ks.begin();
  th_k_iter = th_k.begin();
  for(int k=0; k<K; ++k){

    sum = 0;
    for(int s=0; s<S; ++s){
      // equivalent: beta_w[k*G+g] += th_k_val*h_ks_val;
      sum += (*th_k_iter)*(*h_ks_iter);
      // equivalent: h_ks_val = h_ks[s*K+k];
      h_ks_iter += K; if(s==S-1){h_ks_iter -= K*S;}
    }
    ++h_ks_iter;
    ++th_k_iter;

    for(int g=0; g<G; ++g){
      *beta_w_iter = sum;
      ++beta_w_iter;
    }
  }

  // calculate beta_h
  // equivalent: beta_h_ks[k,]= beta_H_0 + sum(w_gk[,k]*th_k[k] + lamda)
  w_gk_iter = w_gk.begin();
  th_k_iter = th_k.begin();
  for(int k=0; k<K; ++k){
    sum = 0;

    for(int g=0; g<G; ++g){
      sum += (*w_gk_iter)*(*th_k_iter)+lambda;
      // equivalent: w_gk_val = w_gk[k*G+g];
      ++w_gk_iter;
    }
    // equivalent: th_k_val = th_k[k];
    ++th_k_iter;

    for(int s=0; s<S; ++s){
      *beta_h_iter = sum;
      beta_h_iter += K; if(s==S-1){beta_h_iter-=K*S;}
    }
    ++beta_h_iter;
  }


  w_gk_iter = w_gk.begin();
  h_ks_iter = h_ks.begin();
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        // equivalent: beta_th[k] += w_gk_val*h_ks_val;
        *beta_th_iter += (*w_gk_iter)*(*h_ks_iter);
        ++w_gk_iter;
      }
      ++h_ks_iter;
      ++beta_th_iter;
      if(k==K-1){beta_th_iter-=K;}
    }
    w_gk_iter -= G*K;
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector retrofit_step3_alpha_numerator(NumericVector W_gk, NumericVector TH_k, NumericVector H_ks, float lambda) {
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
        phi_a_gks[idx] = (W_gk[(k*G+g)]*TH_k[k] + lambda)*H_ks[(s*K+k)];
        idx++;
      }
    }
  }

  return phi_a_gks;
}

// [[Rcpp::export]]
NumericVector retrofit_step3_alpha_denominator(NumericVector phi_a_gks) {
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

// [[Rcpp::export]]
NumericVector retrofit_step3_alpha(NumericVector W_gk, NumericVector TH_k, NumericVector H_ks, float lambda) {
  /* Equivalent code
   * for(s in 1:S){
   *  for(k in 1:K){
   *    phi_a_gks[,k,s]=((W_gk[,k] * TH_k[k]) +lamda)* H_ks[k,s]
   *  }
   * }
   * for(s in 1:S){
   *  for(v in 1:G){
   *    phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
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
        phi_a_gks[idx] = (W_gk[(k*G+g)]*TH_k[k] + lambda)*H_ks[(s*K+k)];
        idx++;
      }
    }
  }
  
  //Calculate second dimension sums
  idx = 0;
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


// [[Rcpp::export]]
NumericVector retrofit_step3_beta(NumericVector W_gk, NumericVector TH_k, float lambda) {
  /* Equivalent code
   * for(k in 1:K){
   *  for(v in 1:G){
   *    if((W_gk[v,k]*TH_k[k] + lamda)==0){
   *      phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
   *    } else {
   *      phi_b_gk[v,k]=(W_gk[v,k]*TH_k[k])/(W_gk[v,k]*TH_k[k] + lamda)
   *    }
   *  }
   * }
   */
  //dimensions
  // W_gk:(G,K), TH_K:(K)
  NumericVector dim_w = W_gk.attr("dim");
  int G = dim_w[0];
  int K = dim_w[1];
  
  int idx = 0;
  float w = 0;
  float th = 0;
  NumericVector phi_b_gk (G*K);
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      w = W_gk[(k*G+g)];
      th = TH_k[k];
      if((w*th + lambda)==0){
        phi_b_gk[idx]=1; // to avoid numerical error of 0/0 when lambda=0
      } else {
        phi_b_gk[idx]=(w*th)/(w*th + lambda);
      }
      idx++;
    }
  }
  
  return phi_b_gk;
}

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
  NumericVector phi_a_gks (G*K*S);
  NumericVector::iterator phi_a_gks_iter = phi_a_gks.begin();
  NumericVector::iterator W_gk_iter = W_gk.begin();
  NumericVector::iterator H_ks_iter = H_ks.begin();
  NumericVector::iterator TH_k_iter = TH_k.begin();
  // for(int s=0; s<S; ++s){
  //   for(int k=0; k<K; ++k){
  //     for(int g=0; g<G; ++g){
  //       *phi_a_gks_iter = ((*W_gk_iter)*(*TH_k_iter) + lambda)*(*H_ks_iter);
  //       ++phi_a_gks_iter;
  // 
  //       // equivalent: W_gk[k*G+g];
  //       ++W_gk_iter; if(g==G-1 && k==K-1) {W_gk_iter -= K*G;}
  // 
  //       // equivalent: H_ks[s*K+k];
  //       if(g==G-1){++H_ks_iter;}
  // 
  //       // equivalent: TH_k[k];
  //       if(g==G-1) {++TH_k_iter;} if(g==G-1 && k==K-1) {TH_k_iter -= K;}
  //     }
  //   }
  // }

  NumericVector mul_gk (G*K);
  NumericVector::iterator mul_gk_iter = mul_gk.begin();
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      *mul_gk_iter = (*W_gk_iter)*(*TH_k_iter) + lambda;
      ++mul_gk_iter;
      ++W_gk_iter;
    }
    ++TH_k_iter;
  }
  mul_gk_iter = mul_gk.begin();
  float H_ks_val;
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      H_ks_val = *H_ks_iter;
      ++H_ks_iter;
  
      for(int g=0; g<G; ++g){
        *phi_a_gks_iter = (*mul_gk_iter)*H_ks_val;
        ++phi_a_gks_iter;
  
        // equivalent: W_gk[k*G+g];
        ++mul_gk_iter; if(g==G-1 && k==K-1) {mul_gk_iter -= K*G;}
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

  NumericVector sums (G*S);
  NumericVector::iterator sums_iter = sums.begin();
  NumericVector::iterator phi_a_gks_iter = phi_a_gks.begin();
  //Calculate second dimension sums
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        //equivalent: sums[(s*G+g)] += phi_a_gks[idx];
        *sums_iter += *phi_a_gks_iter;
        ++sums_iter; if(g==G-1 && k!=K-1) {sums_iter -= G;}
        ++phi_a_gks_iter;
      }
    }
  }
  
  sums_iter = sums.begin();
  phi_a_gks_iter = phi_a_gks.begin();
  //Calculate denominator part
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        *phi_a_gks_iter = *phi_a_gks_iter/(*sums_iter);
        ++sums_iter; if(g==G-1 && k!=K-1) {sums_iter -= G;}
        ++phi_a_gks_iter;
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
  NumericVector phi_a_gks (G*K*S);
  NumericVector::iterator phi_a_gks_iter = phi_a_gks.begin();
  NumericVector::iterator W_gk_iter = W_gk.begin();
  NumericVector::iterator H_ks_iter = H_ks.begin();
  NumericVector::iterator TH_k_iter = TH_k.begin();
  NumericVector mul_gk (G*K);
  NumericVector::iterator mul_gk_iter = mul_gk.begin();
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      *mul_gk_iter = (*W_gk_iter)*(*TH_k_iter) + lambda;
      ++mul_gk_iter;
      ++W_gk_iter;
    }
    ++TH_k_iter;
  }
  mul_gk_iter = mul_gk.begin();
  float H_ks_val;
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      H_ks_val = *H_ks_iter;
      ++H_ks_iter;
      
      for(int g=0; g<G; ++g){
        *phi_a_gks_iter = (*mul_gk_iter)*H_ks_val;
        ++phi_a_gks_iter;
        
        // equivalent: W_gk[k*G+g];
        ++mul_gk_iter; if(g==G-1 && k==K-1) {mul_gk_iter -= K*G;}
      }
    }
  }
  
  NumericVector sums (G*S);
  NumericVector::iterator sums_iter = sums.begin();
  phi_a_gks_iter = phi_a_gks.begin();
  //Calculate second dimension sums
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        //equivalent: sums[(s*G+g)] += phi_a_gks[idx];
        *sums_iter += *phi_a_gks_iter;
        ++sums_iter; if(g==G-1 && k!=K-1) {sums_iter -= G;}
        ++phi_a_gks_iter;
      }
    }
  }
  
  sums_iter = sums.begin();
  phi_a_gks_iter = phi_a_gks.begin();
  //Calculate denominator part
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      for(int g=0; g<G; ++g){
        *phi_a_gks_iter = *phi_a_gks_iter/(*sums_iter);
        ++sums_iter; if(g==G-1 && k!=K-1) {sums_iter -= G;}
        ++phi_a_gks_iter;
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

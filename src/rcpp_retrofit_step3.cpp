#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void retrofit_decomposition_step3_alpha(List distributions, 
                                        double lambda,
                                        NumericVector dim,
                                        NumericVector &phi_a_gks) {
  /*
   * Equivalent code
   * for(s in 1:S){
   *  for(k in 1:K){
   *    phi_a_gks[,k,s]=((w_gk[,k] * th_k[k]) +lamda)* h_ks[k,s]
   *  }
   * }
   * for(s in 1:S){
   *  for(v in 1:G){
   *    phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
   *  }
   * }
   */
  
  //dimensions
  int G = dim[0];
  int K = dim[1];
  int S = dim[2];
  
  NumericVector w_gk = distributions["w_gk"];
  NumericVector h_ks = distributions["h_ks"];
  NumericVector th_k = distributions["th_k"];
  //Calculate numerator part
  NumericVector::iterator phi_a_gks_iter = phi_a_gks.begin();
  NumericVector::iterator w_gk_iter = w_gk.begin();
  NumericVector::iterator h_ks_iter = h_ks.begin();
  NumericVector::iterator th_k_iter = th_k.begin();
  NumericVector mul_gk (G*K);
  NumericVector::iterator mul_gk_iter = mul_gk.begin();
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      *mul_gk_iter = (*w_gk_iter)*(*th_k_iter) + lambda;
      ++mul_gk_iter;
      ++w_gk_iter;
    }
    ++th_k_iter;
  }
  // numerator variables
  mul_gk_iter = mul_gk.begin();
  double h_ks_val;
  // denominator variables
  NumericVector sums (G*S);
  NumericVector::iterator sums_iter = sums.begin();
  phi_a_gks_iter = phi_a_gks.begin();
  
  for(int s=0; s<S; ++s){
    for(int k=0; k<K; ++k){
      h_ks_val = *h_ks_iter;
      ++h_ks_iter;
      
      for(int g=0; g<G; ++g){
        // numerator part
        *phi_a_gks_iter = (*mul_gk_iter)*h_ks_val;
        // ++phi_a_gks_iter;
        // equivalent: w_gk[k*G+g];
        ++mul_gk_iter; if(g==G-1 && k==K-1) {mul_gk_iter -= K*G;}
        
        // denominator part
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
}


// [[Rcpp::export]]
void retrofit_decomposition_step3_beta(List distributions,
                                       double lambda,
                                       NumericVector dim,
                                       NumericVector &phi_b_gk) {
  /* Equivalent code
   * for(k in 1:K){
   *  for(v in 1:G){
   *    if((w_gk[v,k]*th_k[k] + lamda)==0){
   *      phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
   *    } else {
   *      phi_b_gk[v,k]=(w_gk[v,k]*th_k[k])/(w_gk[v,k]*th_k[k] + lamda)
   *    }
   *  }
   * }
   */
  //dimensions
  int G = dim[0];
  int K = dim[1];
  
  int idx = 0;
  double w = 0;
  double th = 0;
  
  NumericVector w_gk = distributions["w_gk"];
  NumericVector th_k = distributions["th_k"];
  for(int k=0; k<K; ++k){
    for(int g=0; g<G; ++g){
      w = w_gk[(k*G+g)];
      th = th_k[k];
      if((w*th + lambda)==0){
        phi_b_gk[idx]=1; // to avoid numerical error of 0/0 when lambda=0
      } else {
        phi_b_gk[idx]=(w*th)/(w*th + lambda);
      }
      idx++;
    }
  }
}


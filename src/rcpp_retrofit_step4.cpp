#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List retrofit_step4_alpha_calculation(NumericVector x_gs, 
                                     NumericVector phi_a_gks, 
                                     NumericVector phi_b_gk, 
                                     float alpha_w_0,
                                     float alpha_h_0,
                                     float alpha_th_0) {
  
  return List();
}

// [[Rcpp::export]]
List retrofit_step4_beta_calculation(NumericVector W_gk,
                                     NumericVector H_ks,
                                     NumericVector TH_k,
                                     float beta_w_0,
                                     float beta_h_0,
                                     float beta_th_0,
                                     float lambda,
                                     int t)
{
  return List();
}

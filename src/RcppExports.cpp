// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// retrofit_step2_rgamma
NumericVector retrofit_step2_rgamma(NumericVector shapes, NumericVector rates);
RcppExport SEXP _retrofit_retrofit_step2_rgamma(SEXP shapesSEXP, SEXP ratesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shapes(shapesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rates(ratesSEXP);
    rcpp_result_gen = Rcpp::wrap(retrofit_step2_rgamma(shapes, rates));
    return rcpp_result_gen;
END_RCPP
}
// retrofit_step3_alpha_numerator
void retrofit_step3_alpha_numerator(NumericVector W_gk, NumericVector TH_k, NumericVector H_ks, float lambda, NumericVector out_phi_a_gks);
RcppExport SEXP _retrofit_retrofit_step3_alpha_numerator(SEXP W_gkSEXP, SEXP TH_kSEXP, SEXP H_ksSEXP, SEXP lambdaSEXP, SEXP out_phi_a_gksSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type W_gk(W_gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type TH_k(TH_kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type H_ks(H_ksSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type out_phi_a_gks(out_phi_a_gksSEXP);
    retrofit_step3_alpha_numerator(W_gk, TH_k, H_ks, lambda, out_phi_a_gks);
    return R_NilValue;
END_RCPP
}
// retrofit_step3_alpha_denominator
void retrofit_step3_alpha_denominator(NumericVector out_phi_a_gks);
RcppExport SEXP _retrofit_retrofit_step3_alpha_denominator(SEXP out_phi_a_gksSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type out_phi_a_gks(out_phi_a_gksSEXP);
    retrofit_step3_alpha_denominator(out_phi_a_gks);
    return R_NilValue;
END_RCPP
}
// retrofit_step3_alpha
void retrofit_step3_alpha(NumericVector W_gk, NumericVector TH_k, NumericVector H_ks, float lambda, NumericVector out_phi_a_gks);
RcppExport SEXP _retrofit_retrofit_step3_alpha(SEXP W_gkSEXP, SEXP TH_kSEXP, SEXP H_ksSEXP, SEXP lambdaSEXP, SEXP out_phi_a_gksSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type W_gk(W_gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type TH_k(TH_kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type H_ks(H_ksSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type out_phi_a_gks(out_phi_a_gksSEXP);
    retrofit_step3_alpha(W_gk, TH_k, H_ks, lambda, out_phi_a_gks);
    return R_NilValue;
END_RCPP
}
// retrofit_step3_beta
void retrofit_step3_beta(NumericVector W_gk, NumericVector TH_k, float lambda, NumericVector out_phi_b_gk);
RcppExport SEXP _retrofit_retrofit_step3_beta(SEXP W_gkSEXP, SEXP TH_kSEXP, SEXP lambdaSEXP, SEXP out_phi_b_gkSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type W_gk(W_gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type TH_k(TH_kSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type out_phi_b_gk(out_phi_b_gkSEXP);
    retrofit_step3_beta(W_gk, TH_k, lambda, out_phi_b_gk);
    return R_NilValue;
END_RCPP
}
// retrofit_step4_alpha_calculation
List retrofit_step4_alpha_calculation(NumericVector x_gs, NumericVector phi_a_gks, NumericVector phi_b_gk, float alpha_w_0, float alpha_h_0, float alpha_th_0);
RcppExport SEXP _retrofit_retrofit_step4_alpha_calculation(SEXP x_gsSEXP, SEXP phi_a_gksSEXP, SEXP phi_b_gkSEXP, SEXP alpha_w_0SEXP, SEXP alpha_h_0SEXP, SEXP alpha_th_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_gs(x_gsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_a_gks(phi_a_gksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi_b_gk(phi_b_gkSEXP);
    Rcpp::traits::input_parameter< float >::type alpha_w_0(alpha_w_0SEXP);
    Rcpp::traits::input_parameter< float >::type alpha_h_0(alpha_h_0SEXP);
    Rcpp::traits::input_parameter< float >::type alpha_th_0(alpha_th_0SEXP);
    rcpp_result_gen = Rcpp::wrap(retrofit_step4_alpha_calculation(x_gs, phi_a_gks, phi_b_gk, alpha_w_0, alpha_h_0, alpha_th_0));
    return rcpp_result_gen;
END_RCPP
}
// retrofit_step4_beta_calculation
List retrofit_step4_beta_calculation(NumericVector W_gk, NumericVector H_ks, NumericVector TH_k, float beta_w_0, float beta_h_0, float beta_th_0, float lambda);
RcppExport SEXP _retrofit_retrofit_step4_beta_calculation(SEXP W_gkSEXP, SEXP H_ksSEXP, SEXP TH_kSEXP, SEXP beta_w_0SEXP, SEXP beta_h_0SEXP, SEXP beta_th_0SEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type W_gk(W_gkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type H_ks(H_ksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type TH_k(TH_kSEXP);
    Rcpp::traits::input_parameter< float >::type beta_w_0(beta_w_0SEXP);
    Rcpp::traits::input_parameter< float >::type beta_h_0(beta_h_0SEXP);
    Rcpp::traits::input_parameter< float >::type beta_th_0(beta_th_0SEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(retrofit_step4_beta_calculation(W_gk, H_ks, TH_k, beta_w_0, beta_h_0, beta_th_0, lambda));
    return rcpp_result_gen;
END_RCPP
}
// retrofit_step5_parameter_estimation
NumericVector retrofit_step5_parameter_estimation(NumericVector original, NumericVector update, float rho);
RcppExport SEXP _retrofit_retrofit_step5_parameter_estimation(SEXP originalSEXP, SEXP updateSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type original(originalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type update(updateSEXP);
    Rcpp::traits::input_parameter< float >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(retrofit_step5_parameter_estimation(original, update, rho));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_second_dim_sum
NumericVector rcpp_second_dim_sum(NumericVector v);
RcppExport SEXP _retrofit_rcpp_second_dim_sum(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_second_dim_sum(v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_retrofit_retrofit_step2_rgamma", (DL_FUNC) &_retrofit_retrofit_step2_rgamma, 2},
    {"_retrofit_retrofit_step3_alpha_numerator", (DL_FUNC) &_retrofit_retrofit_step3_alpha_numerator, 5},
    {"_retrofit_retrofit_step3_alpha_denominator", (DL_FUNC) &_retrofit_retrofit_step3_alpha_denominator, 1},
    {"_retrofit_retrofit_step3_alpha", (DL_FUNC) &_retrofit_retrofit_step3_alpha, 5},
    {"_retrofit_retrofit_step3_beta", (DL_FUNC) &_retrofit_retrofit_step3_beta, 4},
    {"_retrofit_retrofit_step4_alpha_calculation", (DL_FUNC) &_retrofit_retrofit_step4_alpha_calculation, 6},
    {"_retrofit_retrofit_step4_beta_calculation", (DL_FUNC) &_retrofit_retrofit_step4_beta_calculation, 7},
    {"_retrofit_retrofit_step5_parameter_estimation", (DL_FUNC) &_retrofit_retrofit_step5_parameter_estimation, 3},
    {"_retrofit_rcpp_second_dim_sum", (DL_FUNC) &_retrofit_rcpp_second_dim_sum, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_retrofit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

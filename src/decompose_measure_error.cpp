#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
double decompose_compute_error_two_norm(NumericVector original, NumericVector inferred) {
  double sum = 0;
  NumericVector::iterator original_iter = original.begin();
  NumericVector::iterator inferred_iter = inferred.begin();
  int length = original.length();
  
  for(int i=0; i<length; ++i){
    sum += abs(((*inferred_iter)/(*original_iter)) - 1);
    ++original_iter;
    ++inferred_iter;
  }
  return sum;
}

// [[Rcpp::export]]
double decompose_compute_error_mat_norm(NumericVector original, NumericVector inferred) {
  double max = 0;
  NumericVector::iterator original_iter = original.begin();
  NumericVector::iterator inferred_iter = inferred.begin();
  int length = original.length();
  
  for(int i=0; i<length; ++i){
    if(max < abs(((*inferred_iter)/(*original_iter)) - 1)){
      max = abs(((*inferred_iter)/(*original_iter)) - 1);
    }
    ++original_iter;
    ++inferred_iter;
  }
  return max;
}

// [[Rcpp::export]]
void decompose_update_original(NumericVector original, NumericVector inferred) {
  NumericVector::iterator inferred_iter = inferred.begin();
  NumericVector::iterator original_iter = original.begin();
  int length = inferred.length();
  
  for(int i=0; i<length; ++i){
    *original_iter = *inferred_iter;
    ++inferred_iter;
    ++original_iter;
  }
}

// [[Rcpp::export]]
double decompose_compute_and_update_error_two_norm(NumericVector original, NumericVector inferred) {
  double sum = 0;
  NumericVector::iterator original_iter = original.begin();
  NumericVector::iterator inferred_iter = inferred.begin();
  int length = original.length();
  
  for(int i=0; i<length; ++i){
    sum += abs(((*inferred_iter)/(*original_iter)) - 1);
    *original_iter = *inferred_iter; //update original components as a record for next iteration
    ++original_iter;
    ++inferred_iter;
  }
  return sum;
}

// [[Rcpp::export]]
double decompose_compute_and_update_error_mat_norm(NumericVector original, NumericVector inferred) {
  double max = 0;
  NumericVector::iterator original_iter = original.begin();
  NumericVector::iterator inferred_iter = inferred.begin();
  int length = original.length();
  
  for(int i=0; i<length; ++i){
    if(max < abs(((*inferred_iter)/(*original_iter)) - 1)){
      max = abs(((*inferred_iter)/(*original_iter)) - 1);
    }
    *original_iter = *inferred_iter; //update original components as a record for next iteration
    ++original_iter;
    ++inferred_iter;
  }
  return max;
}
#' RETROFIT 
#' 
#' @description The main algorithm
#'
#' @param x Matrix or Array with two dimensions (GeneExpressions, Spots). 
#'          This is the main spatial transciptomics data.
#' @param ref_w A Matrix or Array with two dimensions (GeneExpressions, Cell types).
#' @param iterations integer: The number of maximum iterations to be executed
#' @param tolerance double: tolerance factor for convergence of the algorithm 
#' @param L integer: The number of components to be decomposed
#' @param K integer: The number of cell types to be selected
#' @param alpha_w_0 double: Variational initial parameter for vector alpha_w
#' @param beta_w_0 double:  Variational initial parameter for vector beta_w
#' @param alpha_h_0 double: Variational initial parameter for vector alpha_h
#' @param beta_h_0 double:  Variational initial parameter for vector beta_h
#' @param alpha_th_0 double:Variational initial parameter for vector alpha_th
#' @param beta_th_0 double: Variational initial parameter for vector beta_th
#' @param lambda double: Background expression profile control
#' @param kappa double: Learning rate factor
#' @param seed double: Random variable seed in case the output should be deterministic
#' @param plot boolean: Plot relative errors
#'
#' @return A list of decomposed vectors that contains
#' \itemize{
#' \item w: 2d array with GeneExpressions, Components
#' \item h: 2d array with Components, Spots
#' \item th: an array with Components
#' \item w_match: 2d array filtered from 'w' with columns labeled as selected cell types (from ref_w)
#' \item h_match: 2d array filtered from 'h' with rows labeled as selected cell types (from ref_w)
#' }
#'
#'@examples
#'x=read.csv(in_path, row.names=1)
#'result = retrofit_decompose(x)
#'@seealso papers reference
#'@export
RetrofitMain <- function(x,
                         ref_w,
                         L           = 16,
                         K           = 8,
                         alpha_w_0   = 0.05, 
                         beta_w_0    = 0.0001, 
                         alpha_h_0   = 0.2,
                         beta_h_0    = 0.2,
                         alpha_th_0  = 1.25,
                         beta_th_0   = 10,
                         lambda      = 0.01,
                         kappa       = 0.5,
                         iterations  = 4000, 
                         tolerance   = 1e-3,
                         seed        = NULL,
                         plot_convergence = FALSE) {
  # Decompose
  ret = RetrofitDecompose(x,
                          L           = L,
                          alpha_w_0   = alpha_w_0, 
                          beta_w_0    = beta_w_0, 
                          alpha_h_0   = alpha_h_0,
                          beta_h_0    = beta_h_0,
                          alpha_th_0  = alpha_th_0,
                          beta_th_0   = beta_th_0,
                          lambda      = lambda,
                          kappa       = kappa,
                          iterations  = iterations, 
                          tolerance   = tolerance,
                          seed        = seed,
                          plot_convergence = plot_convergence)
  w = ret$w
  h = ret$h
  th= ret$th
  # Match
  ret = RetrofitMatch(ref_w, w, h, K)
  w_match = ret$w
  h_match = ret$h
  result <- list(w=w, h=h, th=th, w_match=w_match, h_match=h_match)
  return(result)
}
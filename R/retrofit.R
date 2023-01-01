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
#'iterations = 10
#'L = 16
#'K = 8
#'@seealso papers reference
#'@export
retrofit <- function(x,
                     sc_ref      = NULL,
                     marker_ref   = NULL,
                     L           = 16,
                     K           = 8,
                     iterations  = 4000,
                     lambda      = 0.01,
                     seed        = NULL,
                     alpha_w_0   = 0.05, 
                     beta_w_0    = 0.0001, 
                     alpha_h_0   = 0.2,
                     beta_h_0    = 0.2,
                     alpha_th_0  = 1.25,
                     beta_th_0   = 10,
                     kappa       = 0.5,
                     verbose     = FALSE) {
  # Decompose
  decomposed = decompose(x,
                         L           = L,
                         iterations  = iterations, 
                         lambda      = lambda,
                         seed        = seed,
                         alpha_w_0   = alpha_w_0, 
                         beta_w_0    = beta_w_0, 
                         alpha_h_0   = alpha_h_0,
                         beta_h_0    = beta_h_0,
                         alpha_th_0  = alpha_th_0,
                         beta_th_0   = beta_th_0,
                         kappa       = kappa,
                         verbose     = verbose)
  ret = list(
    decompose = decomposed
  )
  
  # Annotate
  if(!is.null(sc_ref)){
    annotated = annotateWithCorrelations(sc_ref   = sc_ref, 
                                         K        = K, 
                                         decomp_w = decomposed$w, 
                                         decomp_h = decomposed$h)
    ret$annotateWithCorrelations = annotated 
  }
  if(!is.null(marker_ref)){
    annotated = annotateWithMarkers(marker_ref = marker_ref, 
                                    K          = K, 
                                    decomp_w   = decomposed$w, 
                                    decomp_h   = decomposed$h)
    ret$annotateWithMarkers = annotated
  }
  return(ret)
}
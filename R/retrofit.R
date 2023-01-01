#' RETROFIT 
#' 
#' @description The main algorithm
#' @param x           A matrix or array with dimension (GeneExpressions, Spots). This is the main spatial transciptomics data.
#' @param sc_ref      A matrix or array with two dimensions (GeneExpressions, Cell types).
#' @param marker_ref  A list with (keys, values) = (cell types, an array of genes).
#' @param L           integer (default:16)    The number of components to be decomposed
#' @param K integer: The number of cell types to be selected
#' @param iterations  integer (default:4000)  The number of maximum iterations to be executed
#' @param lambda      double  (default:0.01)  Background expression profile control
#' @param seed        double  (default:NULL)  Random variable seed in case the output should be deterministic
#' @param alpha_w_0   double  (default:0.05)  Variational initial parameter for vector alpha_w
#' @param beta_w_0    double  (default:0.0001)Variational initial parameter for vector beta_w
#' @param alpha_h_0   double  (default:0.2)   Variational initial parameter for vector alpha_h
#' @param beta_h_0    double  (default:0.2)   Variational initial parameter for vector beta_h
#' @param alpha_th_0  double  (default:1.25)  Variational initial parameter for vector alpha_th
#' @param beta_th_0   double  (default:10)    Variational initial parameter for vector beta_th
#' @param kappa       double  (default:0.5)   Learning rate factor
#' @param verbose     boolean (default:FALSE)
#'
#' @return A list of decomposed vectors that contains
#' \itemize{
#'    \item decompose: \itemize{ 
#'        \item w:  Decomposed 2d array with GeneExpressions, Components
#'        \item h:  Decomposed 2d array with Components, Spots
#'        \item th: 1d array with Components
#'    }
#'    \item annotateWithCorrelations: \itemize{
#'        \item w:  Filtered 2d array with GeneExpressions, Cell types
#'        \item h:  Filtered2d array with Cell types, Spots
#'    }
#'    \item annotateWithMarkers: \itemize{
#'        \item w:  Filtered 2d array with GeneExpressions, Cell types
#'        \item h:  Filtered2d array with Cell types, Spots
#'    }
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
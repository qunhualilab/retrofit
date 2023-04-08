#' RETROFIT 
#' 
#' @description The main algorithm
#' @param x           A matrix or array with dimension (GeneExpressions, Spots). This is the main spatial transciptomics data.
#' @param sc_ref      A matrix or array with two dimensions (GeneExpressions, Cell types).
#' @param marker_ref  A list with (keys, values) = (cell types, an array of genes).
#' @param L           integer (default:16)    The number of components to be decomposed
#' @param K integer: The number of cell types to be selected
#' @param iterations  integer (default:4000)  The number of maximum iterations to be executed
#' @param init_param  list                    Vatirational initial parameters
#' @param lambda      double  (default:0.01)  Background expression profile control
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
#' utils::data("testSimulationData")
#' iterations  = 10
#' L           = 16
#' K           = 8
#' x           = testSimulationData$extra5_x
#' sc_ref      = testSimulationData$sc_ref
#' 
#' res         = retrofit::retrofit(x, sc_ref=sc_ref, L=L, K=K, iterations=iterations)
#' W           = res$decompose$w
#' W_annotated = res$annotateWithCorrelations$w
#' ranked_cells= res$annotateWithCorrelations$ranked_cells
#'@seealso papers reference
#'@export
retrofit <- function(x,
                     sc_ref      = NULL,
                     marker_ref   = NULL,
                     L           = 16,
                     K           = 8,
                     iterations  = 4000,
                     init_param  = NULL,
                     lambda      = 0.01,
                     kappa       = 0.5,
                     verbose     = FALSE) {
  # Decompose
  decomposed = decompose(x,
                         L           = L,
                         iterations  = iterations, 
                         init_param  = init_param,
                         lambda      = lambda,
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
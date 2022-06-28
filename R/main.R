#' RETROFIT 
#' 
#' @description The main algorithm
#'
#' @param x Matrix(GeneExpressions, Spots): Spatial Transciptomics Data. 
#'
#' @return A list of paths of output data that contain
#' \itemize{
#' \item w: 2d array with GeneExpressions, Components
#' \item h: 2d array with Components, Spots
#' \item th: an array with Components
#' \item match: an array with Components
#' }
#'
#'@examples
#'temporary
#'@seealso papers reference
#'@export
RetrofitMain <- function(x) {
  
  ret = RetrofitDecompose(x, 4000)
  ref_w = NULL
  w = ret$w
  h = ret$h
  RetrofitMatch(ref_w, w, h)
  return(NULL)
}
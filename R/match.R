#' RETROFIT matching algorithm
#' 
#' @description matching algorithm description
#'
#' @param x Matrix(GeneExpressions, Spots): Spatial Transciptomics Data. 
#' @param w
#' @param h
#'
#' @return A list of 
#' \itemize{
#' \item w
#' \item h
#' }
#'
#'@examples
#'@seealso papers reference
#'@export
RetrofitMatch <- function(ref_w, w, h) {
  if(dim(w)[2] != dim(h)[1]){
    stop("dimensions not matched")
  }
  K=dim(ref_w)[2]
  
  correlations = cor(ref_w, w)
  correlations2 = cor(ref_w, w)
  
  ind=rep(NA,K)
  for(i in 1:K){
    r2=which(correlations == max(correlations2), arr.ind=TRUE)[1]
    c2=which(correlations == max(correlations2), arr.ind=TRUE)[2]
    r1=which(correlations2 == max(correlations2), arr.ind=TRUE)[1]
    c1=which(correlations2 == max(correlations2), arr.ind=TRUE)[2]
    ind[r2]=c2
    if(i<K){
      correlations2=correlations2[-r1,-c1]
    }
  }
  w_mod=w[,ind]
  h_mod=h[ind,]
  ret <- list(w=w_mod, h=h_mod, i=ind)
  return(ret)
}
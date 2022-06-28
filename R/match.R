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
  K=dim(w)[2]
  
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
    # match1=which(correlations == max(correlations2), arr.ind=TRUE)
    # match2=which(correlations2 == max(correlations2), arr.ind=TRUE)
    
    # ind[match1[1]]=match1[2]
    # if(i<K){
    #   correlations2=correlations2[-(match2[1]),-(match[2])]
    # }
  }
  w_mod=w[,ind]
  h_mod=h[ind,]
  return <- list(w=w_mod, h=h_mod)
  return(result)
}
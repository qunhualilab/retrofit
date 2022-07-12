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
RetrofitMatch <- function(ref_w, decomp_w, decomp_h, K) {
  if(dim(decomp_w)[2] != dim(decomp_h)[1]){
    stop("dimensions not matched")
  }
  
  # will ref_w always provide cell types?
  cell_types = colnames(ref_w)
  
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w = matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h = matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  
  correlations = cor(ref_w, w)
  correlations2 = cor(ref_w, w)
  print(correlations)
  col_sel=rep(NA,K)
  row_sel=rep(NA,K)
  for(i in 1:K){
    r2=which(correlations == max(correlations2), arr.ind=TRUE)[1]
    c2=which(correlations == max(correlations2), arr.ind=TRUE)[2]
    r1=which(correlations2 == max(correlations2), arr.ind=TRUE)[1]
    c1=which(correlations2 == max(correlations2), arr.ind=TRUE)[2]
    
    print(paste("r2:", r2, cell_types[r2], "r1:", r1, cell_types[r1]))
    
    row_sel[i]=r2
    col_sel[r2]=c2
    
    if(i<K){
      correlations2=correlations2[-r1,-c1]
    }
  }
  
  cell_mod = rep(NA, K)
  for (i in 1:K) {
    cell_mod[i] = cell_types[row_sel[i]]
  }
  
  w_mod = w[,col_sel]
  h_mod = h[col_sel,]
  colnames(w_mod) = cell_mod
  rownames(h_mod) = cell_mod
  
  ret <- list(w=w_mod, h=h_mod, c=cell_mod)
  return(ret)
}
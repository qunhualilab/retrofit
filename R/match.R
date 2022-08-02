#' RETROFIT matching algorithm
#' 
#' @description Match cell types based on correlations with reference. decomp_w   between matching algorithm description
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
RetrofitMatch <- function(ref_w, 
                          decomp_w, 
                          decomp_h, 
                          K) {
  if(dim(decomp_w)[2] != dim(decomp_h)[1]){
    stop("dimensions not matched")
  }
  
  cell_types = colnames(ref_w)
  
  # will ref_w always provide cell types?
  if(is.null(cell_types)){
    col_length = dim(ref_w)[2]
    cell_types = paste('ref_w', array(1:col_length), sep='')
  }
  
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w = matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h = matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  
  # ref_w_normalized = matrix(0, nrow=nrow(ref_w), ncol=ncol(ref_w))
  ref_w_normalized <- ref_w
  ref_w_rowsums = rowSums(ref_w)
  for (i in 1:length(ref_w_rowsums)){
    if(ref_w_rowsums[i] != 0){
      ref_w_normalized[i,] = ref_w[i,]/ref_w_rowsums[i]
    } else {
      ref_w_normalized[i,] = ref_w[i,]
    }
  }
  w_normalized <- w
  w_rowsums = rowSums(w)
  for (i in 1:length(w_rowsums)){
    if(w_rowsums[i] != 0){
      w_normalized[i,] = w[i,]/w_rowsums[i]
    } else {
      w_normalized[i,] = w[i,]
    }
  }
  
  correlations = cor(ref_w_normalized, w_normalized)
  correlations2 = cor(ref_w_normalized, w_normalized)
  
  col_sel=rep(NA,K)
  row_sel=rep(NA,K)
  for(i in 1:K){
    r2=which(correlations == max(correlations2), arr.ind=TRUE)[1]
    c2=which(correlations == max(correlations2), arr.ind=TRUE)[2]
    r1=which(correlations2 == max(correlations2), arr.ind=TRUE)[1]
    c1=which(correlations2 == max(correlations2), arr.ind=TRUE)[2]
    
    row_sel[i]=r2
    col_sel[i]=c2
    
    if(i<K){
      correlations2=correlations2[-r1,-c1]
    }
  }
  # order selections by row numbers
  sel <- data.frame(
    r = row_sel,
    c = col_sel
  )
  sel <- sel[
    with(sel, order(r)),
  ]
  row_sel = sel$r
  col_sel = sel$c
  
  cell_mod = rep(NA, K)
  for (i in 1:K) {
    cell_mod[i] = cell_types[row_sel[i]]
  }
  
  w_mod = w[,col_sel]
  h_mod = h[col_sel,]
  
  print(dim(w_mod))
  colnames(w_mod) = cell_mod
  rownames(h_mod) = cell_mod
  
  ret <- list(w=w_mod, h=h_mod, c=cell_mod)
  return(ret)
}
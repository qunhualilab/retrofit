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
RetrofitMapByCorrelation <- function(ref_cor, 
                                     K,
                                     decomp_w, 
                                     decomp_h) {
  if(dim(decomp_w)[2] != dim(decomp_h)[1]){
    stop("decomp_w and decomp_h dimensions not matched")
  }
  
  cell_types = colnames(ref_cor)
  
  if(is.null(cell_types)){
    col_length = dim(ref_cor)[2]
    cell_types = paste('ref_cor', array(1:col_length), sep='')
  }
  if(length(cell_types)<K){
    warning(paste("cell_types(", length(cell_types), ") are fewer than the mapping target K(", K, "). K is overrided to ", length(cell_types)))
    K = length(cell_types)
  }
  if(dim(decomp_w)[2]<K){
    warning(paste("columns of decomp_w(", dim(decomp_w)[2], ") are fewer than the mapping target K(", K, "). K is overrided to ", dim(decomp_w)[2]))
    K = dim(decomp_w)[2]
  }
  
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w = matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h = matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  rownames(w) = rownames(decomp_w)
  colnames(w) = colnames(decomp_w)
  rownames(h) = rownames(decomp_h)
  colnames(h) = colnames(decomp_h)
  
  # ref_cor_normalized = matrix(0, nrow=nrow(ref_cor), ncol=ncol(ref_cor))
  ref_cor_normalized <- ref_cor
  ref_cor_rowsums = rowSums(ref_cor)
  for (i in 1:length(ref_cor_rowsums)){
    if(ref_cor_rowsums[i] != 0){
      ref_cor_normalized[i,] = ref_cor[i,]/ref_cor_rowsums[i]
    }
  }
  w_normalized <- w
  w_rowsums = rowSums(w)
  
  if(length(w_rowsums) == 0){
    stop("the length of rowsums is 0. the rows of decomposed w may not match with the reference")
  }
  
  for (i in 1:length(w_rowsums)){
    if(w_rowsums[i] != 0){
      w_normalized[i,] = w[i,]/w_rowsums[i]
    }
  }
  
  correlations = cor(ref_cor_normalized, w_normalized)
  correlations2 = cor(ref_cor_normalized, w_normalized)
  
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
  
  colnames(w_mod) = cell_mod
  rownames(h_mod) = cell_mod
  
  ret <- list(w=w_mod, h=h_mod, c=cell_mod)
  return(ret)
}
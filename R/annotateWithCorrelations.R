#' RETROFIT matching algorithm
#' 
#' @description Match cell types based on correlations with reference. decomp_w   between matching algorithm description
#'
#' @param ref_w A Matrix or Array with two dimensions (GeneExpressions, Cell types).
#' @param K integer: The number of cell types to be selected
#' @param decomp_w Matrix(GeneExpressions, Components): Decomposed w matrix
#' @param decomp_h Matrix(Components, Spots): Decomposed h matrix
#'
#' @return A list of 
#' \itemize{
#' \item w
#' \item h
#' }
#'
#'@examples
#'K = 8
#'ref_w=read.csv(paste("../data", "sample_sc_ref.csv", sep="/"), row.names = 1, check.names = FALSE)
#'decomp_w=read.csv(paste("../data/sample_results", "sample_x__decomp_w.csv", sep="/"), row.names = 1, check.names = FALSE)
#'decomp_h=read.csv(paste("../data/sample_results", "sample_x__decomp_h.csv", sep="/"), check.names = FALSE)
#'result = RetrofitMapByCorrelation(sc_ref=ref_w, 
#'                                  K=K,
#'                                  decomp_w = decomp_w,
#'                                  decomp_h = decomp_h)
#'@seealso papers reference
#'@export
annotateWithCorrelations <- function(sc_ref, 
                                     K,
                                     decomp_w, 
                                     decomp_h) {
  testit::assert(!is.null(sc_ref))
  testit::assert(is.numeric(K))
  testit::assert(!is.null(decomp_w))
  testit::assert(!is.null(decomp_h))
  if(dim(decomp_w)[2] != dim(decomp_h)[1]){
    stop("decomp_w and decomp_h dimensions not matched")
  }
  
  cell_types = colnames(sc_ref)
  
  if(is.null(cell_types)){
    col_length = dim(sc_ref)[2]
    cell_types = paste('sc_ref', array(1:col_length), sep='')
  }
  if(length(cell_types)<K){
    warning(paste("cell_types(", length(cell_types), ") are fewer than the mapping target K(", K, "). K is overriden by ", length(cell_types)))
    K = length(cell_types)
  }
  if(dim(decomp_w)[2]<K){
    warning(paste("columns of decomp_w(", dim(decomp_w)[2], ") are fewer than the mapping target K(", K, "). K is overriden by ", dim(decomp_w)[2]))
    K = dim(decomp_w)[2]
  }
  
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w = matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h = matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  rownames(w) = rownames(decomp_w)
  colnames(w) = colnames(decomp_w)
  rownames(h) = rownames(decomp_h)
  colnames(h) = colnames(decomp_h)
  
  # sc_ref_normalized = matrix(0, nrow=nrow(sc_ref), ncol=ncol(sc_ref))
  sc_ref_normalized <- sc_ref
  sc_ref_rowsums = rowSums(sc_ref)
  for (i in 1:length(sc_ref_rowsums)){
    if(sc_ref_rowsums[i] != 0){
      sc_ref_normalized[i,] = sc_ref[i,]/sc_ref_rowsums[i]
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
  
  correlations = cor(sc_ref_normalized, w_normalized)
  correlations2 = cor(sc_ref_normalized, w_normalized)
  
  col_sel=rep(NA,K)
  row_sel=rep(NA,K)
  cor_sel=rep(NA,K) 
  cell_sel=rep(NA, K)
  for(i in 1:K){
    r2=which(correlations == max(correlations2), arr.ind=TRUE)[1]
    c2=which(correlations == max(correlations2), arr.ind=TRUE)[2]
    r1=which(correlations2 == max(correlations2), arr.ind=TRUE)[1]
    c1=which(correlations2 == max(correlations2), arr.ind=TRUE)[2]
    
    row_sel[i]=r2
    col_sel[i]=c2
    cor_sel[i]=correlations[r2,c2]
    cell_sel[i]=cell_types[r2]
    
    if(i<K){
      correlations2=correlations2[-r1,-c1]
    }
  }
  
  # order selections following the reference
  sel <- data.frame(
    r = row_sel,
    c = col_sel
  )
  sel <- sel[
    with(sel, order(r)),
  ]
  row_sel = sel$r
  col_sel = sel$c
  
  w_mod = w[,col_sel]
  h_mod = h[col_sel,]
  cell_mod = rep(NA, K)
  for (i in 1:K) {
    cell_mod[i] = cell_types[row_sel[i]]
  }
  
  colnames(w_mod) = cell_mod
  rownames(h_mod) = cell_mod
  
  ret <- list(w=w_mod, h=h_mod, ranked_cells=cell_sel, ranked_correlations=cor_sel)
  return(ret)
}
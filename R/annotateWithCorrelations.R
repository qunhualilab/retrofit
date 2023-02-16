#' RETROFIT matching algorithm
#' 
#' @description Match cell types based on correlations with reference. decomp_w   between matching algorithm description
#'
#' @param sc_ref A Matrix or Array with two dimensions (GeneExpressions, Cell types).
#' @param K integer: The number of cell types to be selected
#' @param decomp_w Matrix(GeneExpressions, Components): Decomposed w matrix
#' @param decomp_h Matrix(Components, Spots): Decomposed h matrix
#'
#' @return A list of selected components, cells, and correlations
#' \itemize{
#'        \item w:  Filtered 2d array with GeneExpressions, Cell types
#'        \item h:  Filtered2d array with Cell types, Spots
#'        \item ranked_cells: The list of cell names 
#'        \item ranked_correlations: The list of correlations
#' }
#'
#'@examples
#' data("testSimulationData")
#' K            = 10
#' sc_ref       = testSimulationData$sc_ref
#' W            = testSimulationData$decompose$w
#' H            = testSimulationData$decompose$h
#'
#' result       = retrofit::annotateWithCorrelations(sc_ref=sc_ref, K=K, 
#'                                                  decomp_w=W, decomp_h=H)
#' H_annotated  = result$h                                              
#' W_annotated  = result$w
#' ranked_cells = result$ranked_cells
#'@seealso papers reference
#'@export
annotateWithCorrelations <- function(sc_ref, 
                                     K,
                                     decomp_w, 
                                     decomp_h) {
  stopifnot(!is.null(sc_ref))
  stopifnot((is.matrix(sc_ref) || is.array(sc_ref) || is.list(sc_ref)))
  stopifnot(length(dim(sc_ref)) == 2)
  stopifnot(!is.null(rownames(sc_ref)))
  stopifnot(is.numeric(K))
  stopifnot(!is.null(decomp_w))
  stopifnot((is.matrix(decomp_w) || is.array(decomp_w) || is.list(decomp_w)))
  stopifnot(length(dim(decomp_w)) == 2)
  # stopifnot(!is.null(rownames(decomp_w)) && !is.null(colnames(decomp_w)))
  stopifnot(!is.null(decomp_h))
  stopifnot((is.matrix(decomp_h) || is.array(decomp_h) || is.list(decomp_h)))
  stopifnot(length(dim(decomp_h)) == 2)
  # stopifnot(!is.null(rownames(decomp_h)) && !is.null(colnames(decomp_h)))
  stopifnot(dim(decomp_w)[2] == dim(decomp_h)[1])
  stopifnot(dim(decomp_w)[1] == dim(sc_ref)[1])
  
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
  
  correlations = stats::cor(sc_ref_normalized, w_normalized)
  correlations2 = stats::cor(sc_ref_normalized, w_normalized)
  
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
  sel <- data.frame(r=row_sel, c=col_sel)
  sel <- sel[with(sel, order(r)),]
  
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
  
  # weight to proportion
  w_mod_prop <- w_mod
  for(i in 1:nrow(w_mod)){
    w_mod_prop[i,]=w_mod[i,]/sum(w_mod[i,])
  }
  h_mod_prop <- h_mod
  for(i in 1:ncol(h_mod)){
    h_mod_prop[,i]=h_mod[,i]/sum(h_mod[,i])
  }
  
  ret <- list(w=w_mod, 
              h=h_mod, 
              w_prop=w_mod_prop, 
              h_prop=h_mod_prop,
              ranked_cells=cell_sel, 
              ranked_correlations=cor_sel,
              sc_ref_prop=sc_ref_normalized)
  return(ret)
}
#' RETROFIT matching algorithm
#' 
#' @description Match cell types based on correlations with reference. decomp_w   between matching algorithm description
#'
#' @param marker_ref Key-value list: A dictionary of key: cell type, value: GeneExpression list
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
#'decomp_w=read.csv(paste("../data/sample_results", "sample_x__decomp_w.csv", sep="/"), row.names = 1, check.names = FALSE)
#'decomp_h=read.csv(paste("../data/sample_results", "sample_x__decomp_h.csv", sep="/"), check.names = FALSE)
#'marker_ref_d=read.csv(paste("../data", "sample_marker_ref.csv", sep="/"), check.names = FALSE)
#'marker_ref = list()
#'for(r in 1:nrow(marker_ref_d)){
#'  gene = marker_ref_d[[1]][r]
#'  cell_type = marker_ref_d[[2]][r]
#'  if(is.null(marker_ref[[cell_type]])){
#'    marker_ref[[cell_type]] = c()
#'  }
#'  marker_ref[[cell_type]] = c(marker_ref[[cell_type]], gene)
#'}
#'result = RetrofitMapByMarkers(marker_ref = marker_ref, 
#'                              K=K,
#'                              decomp_w = decomp_w,
#'                              decomp_h = decomp_h)
#'@seealso papers reference
#'@export
annotateWithMarkers <- function(marker_ref, 
                                K,
                                decomp_w, 
                                decomp_h) {
  if(dim(decomp_w)[2] != dim(decomp_h)[1]){
    stop("decomp_w and decomp_h dimensions not matched")
  }
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w = matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h = matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  rownames(w) = rownames(decomp_w)
  colnames(w) = colnames(decomp_w)
  rownames(h) = rownames(decomp_h)
  colnames(h) = colnames(decomp_h)
  
  cell_types = names(marker_ref)
  all_genes = unlist(marker_ref)
  
  if(length(cell_types)<K){
    warning(paste("cell_types(", length(cell_types), ") are fewer than the mapping target K(", K, ")."))
    K = length(cell_types)
  }
  if(dim(decomp_w)[2]<K){
    warning(paste("columns of decomp_w(", dim(decomp_w)[2], ") are fewer than the mapping target K(", K, "). K is overriden by ", dim(decomp_w)[2]))
    K = dim(decomp_w)[2]
  }
  
  gene_sums=matrix(NA, nrow=length(cell_types), ncol=ncol(decomp_w))
  rownames(gene_sums)=cell_types
  
  w_normed <- w[rownames(w) %in% all_genes,]
  w_rowsums = rowSums(w_normed)
  
  if(length(w_rowsums) == 0){
    stop("the length of rowsums is 0. the rows of decomposed w may not match with the reference")
  }
  
  for (i in 1:length(w_rowsums)){
    if(w_rowsums[[i]] != 0){
      w_normed[i,] = w_normed[i,]/w_rowsums[i]
    }
  }
  
  for(r in 1:nrow(gene_sums)){
    cell = cell_types[r]
    genes = marker_ref[[cell]]
    genes = unique(genes)
    
    for(c in 1:ncol(gene_sums)){
      w_matched_values = w_normed[rownames(w_normed) %in% genes, c]
      gene_sums[r,c] = sum(w_matched_values)/length(genes)
    }
  }
  
  sums <- gene_sums
  sums2 <- gene_sums
  
  col_sel=rep(NA,K)
  row_sel=rep(NA,K)
  for(i in 1:K){
    r2=which(sums == max(sums2), arr.ind=TRUE)[1]
    c2=which(sums == max(sums2), arr.ind=TRUE)[2]
    r1=which(sums2 == max(sums2), arr.ind=TRUE)[1]
    c1=which(sums2 == max(sums2), arr.ind=TRUE)[2]
    
    row_sel[i]=r2
    col_sel[i]=c2
    
    if(i<K){
      sums2=sums2[-r1,-c1]
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
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
#' utils::data("testSimulationData")
#' K            <- 10
#' marker_ref   <- testSimulationData$marker_ref
#' W            <- testSimulationData$decompose$w
#' H            <- testSimulationData$decompose$h
#'
#' result       <- retrofit::annotateWithMarkers(marker_ref=marker_ref, K=K, 
#'                                             decomp_w=W, decomp_h=H)
#' H_annotated  <- result$h                                              
#' W_annotated  <- result$w
#' ranked_cells <- result$ranked_cells
#'@seealso papers reference
#'@export
annotateWithMarkers <- function(marker_ref, 
                                K,
                                decomp_w, 
                                decomp_h) {
  stopifnot(!is.null(marker_ref))
  stopifnot(is.list(marker_ref))
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
  
  # copy w, h to 'clear' colnames, rownames of w, h respectively.
  w <- matrix(as.numeric(unlist(decomp_w)), nrow=nrow(decomp_w), ncol=ncol(decomp_w))
  h <- matrix(as.numeric(unlist(decomp_h)), nrow=nrow(decomp_h), ncol=ncol(decomp_h))
  rownames(w) <- rownames(decomp_w)
  colnames(w) <- colnames(decomp_w)
  rownames(h) <- rownames(decomp_h)
  colnames(h) <- colnames(decomp_h)
  
  cell_types <- names(marker_ref)
  all_genes <- unlist(marker_ref)
  
  if(length(cell_types)<K){
    warning(paste("cell_types(", length(cell_types), ") are fewer than the mapping target K(", K, ")."))
    K <- length(cell_types)
  }
  if(dim(decomp_w)[2]<K){
    warning(paste("columns of decomp_w(", dim(decomp_w)[2], ") are fewer than the mapping target K(", K, "). K is overriden by ", dim(decomp_w)[2]))
    K <- dim(decomp_w)[2]
  }
  
  gene_sums <- matrix(NA, nrow=length(cell_types), ncol=ncol(decomp_w))
  rownames(gene_sums) <- cell_types
  
  w_normed <- w[rownames(w) %in% all_genes,]
  w_rowsums <- rowSums(w_normed)
  
  if(length(w_rowsums) == 0){
    stop("the length of rowsums is 0. the rows of decomposed w may not match with the reference")
  }
  
  for (i in seq_along(w_rowsums)){
    if(w_rowsums[[i]] != 0){
      w_normed[i,] <- w_normed[i,]/w_rowsums[i]
    }
  }
  
  for(r in seq_len(nrow(gene_sums))){
    cell <- cell_types[r]
    genes <- marker_ref[[cell]]
    genes <- unique(genes)
    
    for(c in seq_len(ncol(gene_sums))){
      w_matched_values <- w_normed[rownames(w_normed) %in% genes, c]
      gene_sums[r,c] <- sum(w_matched_values)/length(genes)
    }
  }
  
  sums <- gene_sums
  sums2 <- gene_sums
  
  col_sel <- rep(NA,K)
  row_sel <- rep(NA,K)
  cell_sel <- rep(NA, K)
  for(i in seq_len(K)){
    r2 <- which(sums == max(sums2), arr.ind=TRUE)[1]
    c2 <- which(sums == max(sums2), arr.ind=TRUE)[2]
    r1 <- which(sums2 == max(sums2), arr.ind=TRUE)[1]
    c1 <- which(sums2 == max(sums2), arr.ind=TRUE)[2]
    
    row_sel[i] <- r2
    col_sel[i] <- c2
    cell_sel[i] <- cell_types[r2]
    
    if(i<K){
      sums2 <- sums2[-r1,-c1]
    }
  }
  
  # order selections by row numbers
  sel <- data.frame(r=row_sel,c=col_sel)
  sel <- sel[with(sel, order(r)),]
  row_sel <- sel$r
  col_sel <- sel$c
  
  cell_mod <- rep(NA, K)
  for (i in seq_len(K)) {
    cell_mod[i] <- cell_types[row_sel[i]]
  }
  
  w_mod <- w[,col_sel]
  h_mod <- h[col_sel,]
  colnames(w_mod) <- cell_mod
  rownames(h_mod) <- cell_mod
  
  # weight to proportion
  w_mod_prop <- w_mod
  for(i in seq_len(nrow(w_mod))){
    w_mod_prop[i,] <- w_mod[i,]/sum(w_mod[i,])
  }
  h_mod_prop <- h_mod
  for(i in seq_len(ncol(h_mod))){
    h_mod_prop[,i] <- h_mod[,i]/sum(h_mod[,i])
  }
  
  ret <- list(w=w_mod, 
              h=h_mod, 
              w_prop=w_mod_prop, 
              h_prop=h_mod_prop,
              ranked_cells=cell_sel, 
              gene_sums=sums)
  return(ret)
}
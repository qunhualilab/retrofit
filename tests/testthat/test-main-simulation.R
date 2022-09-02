test_that("main-simulation-works", {
  iterations = 10
  file_source = "sample_x"
  ref_file_source = "sample_ref_cor"
  ref_marker_file_source = "sample_ref_marker"
  L = 16
  K = 8
  seed = 1
  
  in_dir = "../../data"
  result_dir = "../../data/sample_results"
  ref_w_path = paste(in_dir, paste(ref_file_source, ".csv", sep=""), sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1, check.names = FALSE)
  ref_marker_path = paste(in_dir, paste(ref_marker_file_source, ".csv", sep=""), sep="/")
  ref_marker_d=read.csv(ref_marker_path, check.names = FALSE)
  ref_marker = list()
  for(r in 1:nrow(ref_marker_d)){
    gene = ref_marker_d[[1]][r]
    cell_type = ref_marker_d[[2]][r]
    if(is.null(ref_marker[[cell_type]])){
      ref_marker[[cell_type]] = c()
    }
    ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  }
  
  in_path = paste(in_dir, paste(file_source, ".csv", sep=""), sep="/")
  x=read.csv(in_path, row.names = 1, check.names = FALSE)
  
  result = Retrofit(x, 
                    ref_cor=ref_w, 
                    ref_marker=ref_marker, 
                    iterations=iterations, 
                    L=L, 
                    K=K,
                    seed=seed)
  
  decomp_h_path       = paste(result_dir, paste(file_source,"__decomp_h.csv", sep=""), sep="/")
  decomp_w_path       = paste(result_dir, paste(file_source,"__decomp_w.csv", sep=""), sep="/")
  decomp_th_path      = paste(result_dir, paste(file_source,"__decomp_th.csv", sep=""), sep="/")
  match_cor_w_path    = paste(result_dir, paste(file_source,"__map_cor_w.csv", sep=""), sep="/")
  match_cor_h_path    = paste(result_dir, paste(file_source,"__map_cor_h.csv", sep=""), sep="/")
  match_marker_w_path = paste(result_dir, paste(file_source,"__map_marker_w.csv", sep=""), sep="/")
  match_marker_h_path = paste(result_dir, paste(file_source,"__map_marker_h.csv", sep=""), sep="/")
  
  decomp_h        = read.csv(decomp_h_path, row.names = 1, check.names = FALSE)
  decomp_w        = read.csv(decomp_w_path, row.names = 1, check.names = FALSE)
  decomp_th       = read.csv(decomp_th_path, row.names = 1, check.names = FALSE)
  match_cor_h     = read.csv(match_cor_h_path, row.names = 1, check.names = FALSE)
  match_cor_w     = read.csv(match_cor_w_path, row.names = 1, check.names = FALSE)
  match_marker_h  = read.csv(match_marker_h_path, row.names = 1, check.names = FALSE)
  match_marker_w  = read.csv(match_marker_w_path, row.names = 1, check.names = FALSE)
  
  expect_true(all.equal(matrix(unlist(decomp_h)), result$h, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(decomp_w)), result$w, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(decomp_th)), result$t, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(match_cor_h)), result$h_cor_map, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(match_cor_w)), result$w_cor_map, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(match_marker_h)), result$h_marker_map, check.attributes = FALSE))
  expect_true(all.equal(matrix(unlist(match_marker_w)), result$w_marker_map, check.attributes = FALSE))
})




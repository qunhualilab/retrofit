RetrofitMainSimulation <- function(folder) {
  setwd("~/Research/retrofit/retrofit/R")
  
  # file_source = "N=20,M=5_loc_X"
  # file_source = "N=10,M=3_loc_X"
  file_source = "extra5_loc_X"
  ref_file_source = "Cerebellum_W_K=10"
  L = 10
  K = 5
  in_dir = "../results"
  out_dir = "../results/local"
  iterations=4000
  
  if(!is.null(folder)){
    out_dir = paste(out_dir, folder, sep="/")
    dir.create(out_dir)
  }
  ref_w_path = paste(in_dir, paste(ref_file_source, ".csv", sep=""), sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  in_path = paste(in_dir, paste(file_source, ".csv", sep=""), sep="/")
  x=read.csv(in_path, row.names = 1)
  
  result = RetrofitMain(x, ref_w, iterations=iterations, L=L, K=K)
  
  decomp_h_path = paste(out_dir, paste(file_source,"__decomp_h.csv", sep=""), sep="/")
  decomp_w_path = paste(out_dir, paste(file_source,"__decomp_w.csv", sep=""), sep="/")
  decomp_th_path = paste(out_dir, paste(file_source,"__decomp_th.csv", sep=""), sep="/")
  match_w_path = paste(out_dir, paste(file_source,"__match_w.csv", sep=""), sep="/")
  match_h_path = paste(out_dir, paste(file_source,"__match_h.csv", sep=""), sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_match, match_h_path)
  write.csv(result$w_match, match_w_path)
  
  paths <- list(decomp_w_path, decomp_h_path, decomp_th_path, match_w_path, match_h_path)
  print(paths)
  return(paths)
}

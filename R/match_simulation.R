MatchSimulation <- function(dir, file) {
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  out_h_file = paste(file,"_out_H.csv", sep="")
  out_h_path = paste(dir, out_h_file, sep="/")
  out_w_file = paste(file,"_out_W.csv", sep="")
  out_w_path = paste(dir, out_w_file, sep="/")
  out_t_file = paste(file,"_out_T.csv", sep="")
  out_t_path = paste(dir, out_t_file, sep="/")
  
  X=read.csv(in_path)
  rownames(X)=X[,1]
  X=as.matrix(X[,-1])
  
  result = RetrofitMatch(X, iterations)
  
  write.csv(result["h"], out_h_path)
  write.csv(result["w"], out_w_path)
  write.csv(result["t"], out_t_path)
  
  print("simulation finished")
  paths <- list(in_path, out_w_path, out_h_path, out_t_path)
  names(paths) <- c("in_path", "out_w_path", "out_h_path", "out_t_path")
  return(paths)
}

RetrofitMatchSimulationLocal <- function() {
  print(paste("working directory: ", getwd()))
  setwd("~/Research/retrofit/retrofit/R")
  
  # in_file = c("N=20,M=5_loc_X", "N=10,M=3_loc_X", "extra5_loc_X")
  in_file = "N=10,M=3_loc_X"
  in_dir = "../results"
  K = 10
  
  ref_w_file = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  decomp_dir = in_dir
  match_dir = in_dir
  decomp_h_file = paste(in_file,"_decomp_H_curved.csv", sep="")
  decomp_h_path = paste(decomp_dir, decomp_h_file, sep="/")
  decomp_w_file = paste(in_file,"_decomp_W_curved.csv", sep="")
  decomp_w_path = paste(decomp_dir, decomp_w_file, sep="/")
  
  match_h_file = paste(in_file,"_match_H_curved.csv", sep="")
  match_h_path = paste(match_dir, match_h_file, sep="/")
  match_w_file = paste(in_file,"_match_W_curved.csv", sep="")
  match_w_path = paste(match_dir, match_w_file, sep="/")
  
  
  w=read.csv(decomp_w_path, row.names = 1)
  h=read.csv(decomp_h_path, row.names = 1)
  
  ret = RetrofitMatchWithRef(ref_w, w, h, K)
  
  write.csv(ret$w, match_w_path)
  write.csv(ret$h, match_h_path)
}

RetrofitMatchSimulationTemp <- function() {
  print(paste("working directory: ", getwd()))
  
  in_r_dir = "~/Research/retrofit/retrofit/results/local/plot"
  in_w_file = "N=20,M=5_loc_W_hat_L=20"
  in_h_file = "N=20,M=5_loc_H_hat_L=20"
  K = 10
  
  w_path = paste(in_r_dir,"/", in_w_file, '.csv', sep="")
  h_path = paste(in_r_dir,"/", in_h_file, '.csv', sep="")
  
  w=read.csv(w_path, row.names = 1)
  h=read.csv(h_path, row.names = 1)
  
  setwd("~/Research/retrofit/retrofit/R")
  in_dir = "../results"
  ref_w_file = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  ret = RetrofitMatchWithRef(ref_w, w, h, K)
  
  write.csv(ret$w, paste(in_dir, "/", "local", "/", in_w_file, '_match.csv', sep=""))
  write.csv(ret$h, paste(in_dir, "/", "local", "/", in_h_file, '_match.csv', sep=""))
}

RetrofitMatchMarkerTemp <- function() {
  print(paste("working directory: ", getwd()))
  
  in_r_dir = "~/Research/retrofit/Simulation backup"
  in_w_file = "A4_X_decomposed_W"
  in_h_file = "A4_X_decomposed_H"
  
  w_path = paste(in_r_dir,"/", in_w_file, '.csv', sep="")
  h_path = paste(in_r_dir,"/", in_h_file, '.csv', sep="")
  
  w=read.csv(w_path, row.names = 1)
  h=read.csv(h_path, row.names = 1)
  
  setwd("~/Research/retrofit/retrofit/R")
  in_dir = "../results"
  ref_w_file = "Colon_Markers.csv"
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  ref_marker = list()
  for(r in 1:nrow(ref_w)){
    gene = ref_w[[1]][r]
    cell_type = ref_w[[2]][r]
    if(is.null(ref_marker[[cell_type]])){
      ref_marker[[cell_type]] = c()
    }
    ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  }
  
  ret = RetrofitMapByMarkers(ref_marker, w, h)
  
  write.csv(ret$w, paste(in_dir, "/", "local", "/", in_w_file, '_match_markers.csv', sep=""))
  write.csv(ret$h, paste(in_dir, "/", "local", "/", in_h_file, '_match_markers.csv', sep=""))
}
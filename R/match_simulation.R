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
  
  result = retrofit(X, iterations)
  
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
  
  in_file = "A3_1554"
  in_dir = "../results"
  decomp_dir = in_dir
  match_dir = in_dir
  decomp_h_file = paste(in_file,"_decomp_H.csv", sep="")
  decomp_h_path = paste(decomp_dir, decomp_h_file, sep="/")
  decomp_w_file = paste(in_file,"_decomp_W.csv", sep="")
  decomp_w_path = paste(decomp_dir, decomp_w_file, sep="/")
  
  match_h_file = paste(in_file,"_match_H.csv", sep="")
  match_h_path = paste(match_dir, match_h_file, sep="/")
  match_w_file = paste(in_file,"_match_W.csv", sep="")
  match_w_path = paste(match_dir, match_w_file, sep="/")
  
  ref_w_file = "A3_1554_ref_W"
  ref_w_file = paste(ref_w_file,".csv", sep="")
  ref_w_path = paste(match_dir, ref_w_file, sep="/")
  
  w=read.csv(decomp_w_path)
  w=w[,-1]
  h=read.csv(decomp_h_path)
  h=h[,-1]
  ref_w=read.csv(ref_w_path)
  ref_w=ref_w[,-1]
  
  print(dim(ref_w))
  print(dim(w))
  print(dim(h))
  
  write.csv(cor(ref_w, w), match_w_path)
  ret = RetrofitMatch(ref_w,w,h)
  write.csv(ret$w, match_w_path)
  write.csv(ret$h, match_h_path)
}
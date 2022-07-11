RetrofitDecomposeSimulation <- function(dir, file, iterations=2) {
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  decomp_h_file = paste(file,"_decomp_H.csv", sep="")
  decomp_h_path = paste(dir, decomp_h_file, sep="/")
  decomp_w_file = paste(file,"_decomp_W.csv", sep="")
  decomp_w_path = paste(dir, decomp_w_file, sep="/")
  decomp_t_file = paste(file,"_decomp_T.csv", sep="")
  decomp_t_path = paste(dir, decomp_t_file, sep="/")
  
  x=read.csv(in_path)
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  result = RetrofitDecompose(x, iterations=iterations)
  
  write.csv(result["h"], decomp_h_path)
  write.csv(result["w"], decomp_w_path)
  write.csv(result["t"], decomp_t_path)
  
  print("simulation finished")
  paths <- list(in_path, decomp_w_path, decomp_h_path, decomp_t_path)
  names(paths) <- c("decomp_in", "decomp_w", "decomp_h", "decomp_t")
  return(paths)
}

RetrofitDecomposeSimulationLocal <- function() {
  # file_source = "N=20,M=5_loc_X"
  file_source = "A3_1554"
  in_dir = "results"
  decomp_dir = in_dir
  iterations=4000
  in_file = paste(file_source,".csv", sep="")
  in_path = paste(in_dir, in_file, sep="/")
  decomp_h_file = paste(file_source,"_decomp_H.csv", sep="")
  decomp_h_path = paste(decomp_dir, decomp_h_file, sep="/")
  decomp_w_file = paste(file_source,"_decomp_W.csv", sep="")
  decomp_w_path = paste(decomp_dir, decomp_w_file, sep="/")
  decomp_t_file = paste(file_source,"_decomp_T.csv", sep="")
  decomp_t_path = paste(decomp_dir, decomp_t_file, sep="/")

  print(in_path)
  x=read.csv(in_path)
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])

  result = RetrofitDecompose(x, iterations=iterations)

  write.csv(result["h"], decomp_h_path)
  write.csv(result["w"], decomp_w_path)
  write.csv(result["t"], decomp_t_path)
  
  print(paste("simulation finished saved at: ", decomp_w_path, decomp_h_path, decomp_t_path, sep = " "))
}
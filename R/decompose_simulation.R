retrofit_simulation <- function(dir, file, iterations=2) {
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

retrofit_simulation_local <- function() {
  in_file = "a3_1554.csv"
  in_dir = "../../results"
  out_dir = "../../results/local"
  out_file = "temp"
  iterations=2
  in_path = paste(in_dir, in_file, sep="/")
  out_h_file = paste(out_file,"_out_H.csv", sep="")
  out_h_path = paste(out_dir, out_h_file, sep="/")
  out_w_file = paste(out_file,"_out_W.csv", sep="")
  out_w_path = paste(out_dir, out_w_file, sep="/")
  out_t_file = paste(out_file,"_out_T.csv", sep="")
  out_t_path = paste(out_dir, out_t_file, sep="/")
  
  X=read.csv(in_path)
  
  result = retrofit(X, iterations)
  
  write.csv(result["h"], out_h_path)
  write.csv(result["w"], out_w_path)
  write.csv(result["t"], out_t_path)
  
  print("simulation finished")
}
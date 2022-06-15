retrofit_simulation_original <- function(dir, file, iterations=2) {
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  out_h_file = paste(file,"_H","_original.csv", sep="")
  out_h_path = paste(dir, out_h_file, sep="/")
  out_w_file = paste(file,"_W","_original.csv", sep="")
  out_w_path = paste(dir, out_w_file, sep="/")
  out_t_file = paste(file,"_T","_original.csv", sep="")
  out_t_path = paste(dir, out_t_file, sep="/")

  X=read.csv(in_path)
  
  result = retrofit_original(X, iterations)

  write.csv(result["h"], out_h_path)
  write.csv(result["w"], out_w_path)
  write.csv(result["t"], out_t_path)

  print("simulation finished")
  paths <- list(in_path, out_w_path, out_h_path, out_t_path)
  names(paths) <- c("in_path", "out_w_path", "out_h_path", "out_t_path")
  return(paths)
}
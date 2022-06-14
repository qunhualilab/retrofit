# rm(list=ls())

a3_1554 <- function() {
  dir_path = "~/Research/retrofit/codes"
  in_origin_path = "A3_counts_G=1554.csv"
  out_h_path = "A3_H_hat_G=1554.csv"
  out_w_path = "A3_W_hat_G=1554.csv"
  out_t_path = "A3_Theta_hat_G=1554.csv"

  X=read.csv(paste(dir_path, in_origin_path, sep="/"))
  ######## X = matrix(X, nrow=1081, ncol=1555)

  result = retrofit(x=X, iterations=400)

  write.csv(result["h"], paste(dir_path, out_h_path, sep="/"))
  write.csv(result["w"], paste(dir_path, out_w_path, sep="/"))
  write.csv(result["t"], paste(dir_path, out_t_path, sep="/"))

  print("simulation finished")
}
test_that("step4_alpha", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  x=matrix(runif(G*S, 0, 10),nrow=G, ncol=S)
  phi_a_gks=array(runif(G*K*S, 0, 10), c(G,K,S))
  phi_b_gk=array(runif(G*K, 0, 10), c(G,K))
  alpha_w_0 = runif(1, 0, 1)
  alpha_h_0 = runif(1, 0, 1)
  alpha_th_0 = runif(1, 0, 1)
  
  # r code
  from = Sys.time()
  alpha_w_gk = matrix(rep(0, G*K) , nrow=G, ncol=K)
  alpha_h_ks = matrix(rep(0, K*S) , nrow=K, ncol=S)
  alpha_th_k = array(rep(0, K) , c(K))
  for(k in 1:K){
    alpha_w_gk[,k] = alpha_w_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k]
    alpha_h_ks[k,] = alpha_h_0 + colSums(x*phi_a_gks[,k,])
    alpha_th_k[k]  = alpha_th_0+ sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
  }
  print(paste('r-code: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  # rcpp code
  from = Sys.time()
  ret = retrofit_step4_alpha_calculation(x, 
                                         phi_a_gks, 
                                         phi_b_gk, 
                                         alpha_w_0,
                                         alpha_h_0,
                                         alpha_th_0)
  alpha_w_gk_new = matrix(ret$w, nrow=G, ncol=K)
  alpha_h_ks_new = matrix(ret$h, nrow=K, ncol=S)
  alpha_th_k_new = array(ret$t , c(K))
  print(paste('rcpp-code: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  
  expect_true(all.equal(alpha_w_gk, alpha_w_gk_new, tolerance=1e-4))
  expect_true(all.equal(alpha_h_ks, alpha_h_ks_new, tolerance=1e-4))
  expect_true(all.equal(alpha_th_k, alpha_th_k_new, tolerance=1e-4))
})

test_that("step4_beta", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  beta_w_0 = runif(1, 0, 1)
  beta_h_0 = runif(1, 0, 1)
  beta_th_0 = runif(1, 0, 1)

  # r code
  from = Sys.time()
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  H_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  TH_k=runif(K, 0, 10)
  lambda = runif(1, 0, 1)
  beta_w_gk=matrix(runif(G*K,0,0.005)+beta_w_0 , nrow=G, ncol=K)
  beta_h_ks=matrix(runif(K*S,0,0.5)+beta_h_0 , nrow=K, ncol=S)
  beta_th_k=runif(K,0,1)+beta_th_0
  for(k in 1:K){
    beta_w_gk[,k]= beta_w_0 + sum(H_ks[k,]*TH_k[k])
    beta_h_ks[k,]= beta_h_0 + sum(W_gk[,k]*TH_k[k] + lambda)
    beta_th_k[k]= beta_th_0 + sum(as.matrix(W_gk[,k]) %*% t(as.matrix(H_ks[k,])))
  }
  print(paste('r-code: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))

  # rcpp code
  from = Sys.time()
  ret = retrofit_step4_beta_calculation(W_gk,
                                        H_ks,
                                        TH_k,
                                        beta_w_0,
                                        beta_h_0,
                                        beta_th_0,
                                        lambda)
  beta_w_gk_new = matrix(ret$w, nrow=G, ncol=K)
  beta_h_ks_new = matrix(ret$h, nrow=K, ncol=S)
  beta_th_k_new = ret$t
  print(paste('rcpp-code: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))

  expect_true(all.equal(beta_w_gk, beta_w_gk_new, tolerance=1e-4))
  expect_true(all.equal(beta_h_ks, beta_h_ks_new, tolerance=1e-4))
  expect_true(all.equal(beta_th_k, beta_th_k_new, tolerance=1e-4))
})
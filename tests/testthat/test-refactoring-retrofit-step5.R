test_that("step5-alpha-updates", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  W_gk_updated=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  H_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  H_ks_updated=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  TH_k=array(runif(K, 0, 10))
  TH_k_updated=array(runif(K, 0, 10))
  rho = runif(1, 0, 1)
  
  alpha_W_gk = array(rep(0,G*K), c(G,K))
  alpha_H_ks = array(rep(0,K*S), c(K,S))
  alpha_TH_k = array(rep(0,K), c(K))
  # r code
  for(k in 1:K){
    alpha_W_gk[,k]= (1-rho)*W_gk[,k] + rho*W_gk_updated[,k]
    alpha_H_ks[k,]= (1-rho)*H_ks[k,] + rho*H_ks_updated[k,]
    alpha_TH_k[k]= (1-rho)*TH_k[k] + rho*TH_k_updated[k]
  }
  
  # rcpp code
  alpha_W_gk_new = array(retrofit_step5_parameter_estimation(W_gk, W_gk_updated, rho), c(G,K))
  alpha_H_ks_new = array(retrofit_step5_parameter_estimation(H_ks, H_ks_updated, rho), c(K,S))
  alpha_TH_k_new = array(retrofit_step5_parameter_estimation(TH_k, TH_k_updated, rho), c(K))
  
  expect_true(all.equal(alpha_W_gk, alpha_W_gk_new))
  expect_true(all.equal(alpha_H_ks, alpha_H_ks_new))
  expect_true(all.equal(alpha_TH_k, alpha_TH_k_new))
})

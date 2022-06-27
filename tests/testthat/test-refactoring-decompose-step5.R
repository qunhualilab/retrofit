test_that("step5-alpha-updates", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  w_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  w_gk_updated=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  h_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  h_ks_updated=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  th_k=array(runif(K, 0, 10))
  th_k_updated=array(runif(K, 0, 10))
  rho = runif(1, 0, 1)
  
  alpha_w_gk = array(rep(0,G*K), c(G,K))
  alpha_h_ks = array(rep(0,K*S), c(K,S))
  alpha_th_k = array(rep(0,K), c(K))
  # r code
  for(k in 1:K){
    alpha_w_gk[,k]= (1-rho)*w_gk[,k] + rho*w_gk_updated[,k]
    alpha_h_ks[k,]= (1-rho)*h_ks[k,] + rho*h_ks_updated[k,]
    alpha_th_k[k]= (1-rho)*th_k[k] + rho*th_k_updated[k]
  }
  
  # rcpp code
  alpha_w_gk_new <- w_gk
  alpha_h_ks_new <- h_ks
  alpha_th_k_new <- th_k
  decompose_step5(alpha_w_gk_new, w_gk_updated, rho)
  decompose_step5(alpha_h_ks_new, h_ks_updated, rho)
  decompose_step5(alpha_th_k_new, th_k_updated, rho)
  
  expect_true(all.equal(alpha_w_gk, alpha_w_gk_new))
  expect_true(all.equal(alpha_h_ks, alpha_h_ks_new))
  expect_true(all.equal(alpha_th_k, alpha_th_k_new))
})

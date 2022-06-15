test_that("step4_alpha_w", {
  G = 4
  K = 3
  S = 2
  x_gs=matrix(runif(G*S, 0, 10),nrow=G, ncol=S)
  phi_a_gks=array(runif(G*K*S, 0, 10), c(G,K,S))
  phi_b_gk=array(runif(G*K, 0, 10), c(G,K))
  alpha_w_0 = runif(1, 0, 1)
  
  alpha_w_gk = matrix(runif(G*K,0,0.5)+alpha_w_0 , nrow=G, ncol=K)
  alpha_w_gk_new <- alpha_w_gk
  
  # r code
  for(k in 1:K){
    alpha_w_gk[,k] = alpha_w_0 + rowSums(x_gs*phi_a_gks[,k,])*phi_b_gk[,k]
  }
  
  # rcpp code
  alpha_w_gk_new = matrix(retrofit_step4_alpha_w(x_gs, phi_a_gks, phi_b_gk, alpha_w_0), nrow=G, ncol=K);
  
  expect_true(all.equal(alpha_w_gk, alpha_w_gk_new, tolerance=1e-4))
})

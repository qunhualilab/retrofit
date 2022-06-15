test_that("step3-alpha-numerator", {
  G = 4
  K = 3
  S = 2
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  H_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  TH_k=array(runif(K, 0, 10))
  lambda = runif(1, 0, 1)
  
  # r code
  phi_a_gks=array(rep(0,G*K*S), c(G,K,S))
  for(s in 1:S){
    for(k in 1:K){
      phi_a_gks[,k,s] = ((W_gk[,k] * TH_k[k]) +lambda)* H_ks[k,s]
    }
  }
  
  # rcpp code
  phi_a_gks_new = array(retrofit_step3_alpha_numerator(W_gk, TH_k, H_ks, lambda), c(G,K,S));
  
  expect_true(all.equal(phi_a_gks, phi_a_gks_new, tolerance=1e-4))
})

test_that("step3-alpha-denominator", {
  G = 4
  K = 3
  S = 2
  sample=array(runif(G*K*S, 0, 10), c(G,K,S))
  
  # r code
  phi_a_gks <- sample
  for(s in 1:S){
    for(v in 1:G){
      phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
    }
  }
  
  # rcpp code
  phi_a_gks_new <- sample
  phi_a_gks_new = array(retrofit_step3_alpha_denominator(phi_a_gks_new), c(G,K,S));
  
  # test: all elements are same
  expect_true(all.equal(phi_a_gks, phi_a_gks_new, tolerance=1e-4))
})

test_that("step3-beta", {
  G = 4
  K = 3
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  TH_k=array(runif(K, 0, 10))
  lambda = runif(1, 0, 1)
  
  # r code
  phi_b_gk=array(rep(0,G*K), c(G,K))
  for(k in 1:K){
    for(v in 1:G){
      if((W_gk[v,k]*TH_k[k] + lambda)==0){
        phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lambda=0
      } else {
        phi_b_gk[v,k]=(W_gk[v,k]*TH_k[k])/(W_gk[v,k]*TH_k[k] + lambda)
      }
    }
  }
  
  # rcpp code
  phi_b_gk_new = array(retrofit_step3_beta(W_gk, TH_k, lambda), c(G,K));
  
  expect_true(all.equal(phi_b_gk, phi_b_gk_new, tolerance=1e-4))
})

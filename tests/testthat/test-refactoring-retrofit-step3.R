test_that("step3-alpha", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  H_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  TH_k=array(runif(K, 0, 10))
  lambda = runif(1, 0, 1)
  phi_a_gks = array(rep(0,G*K*S), c(G,K,S))
  phi_a_gks_new <- phi_a_gks
  
  # r code
  from = Sys.time()
  for(s in 1:S){
    for(k in 1:K){
      phi_a_gks[,k,s] = ((W_gk[,k] * TH_k[k]) +lambda)* H_ks[k,s]
    }
    for(v in 1:G){
      phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
    }
  }
  print(paste('r: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  
  # rcpp code
  from = Sys.time()
  retrofit_step3_alpha(W_gk, TH_k, H_ks, lambda, phi_a_gks_new, c(G,K,S));
  print(paste('rcpp: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  expect_true(all.equal(phi_a_gks, phi_a_gks_new))
})

test_that("step3-beta", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  W_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  TH_k=array(runif(K, 0, 10))
  lambda = runif(1, 0, 1)
  phi_b_gk=array(rep(0,G*K), c(G,K))
  phi_b_gk_new <- phi_b_gk
  
  # r code
  # from = Sys.time()
  for(k in 1:K){
    for(v in 1:G){
      if((W_gk[v,k]*TH_k[k] + lambda)==0){
        phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lambda=0
      } else {
        phi_b_gk[v,k]=(W_gk[v,k]*TH_k[k])/(W_gk[v,k]*TH_k[k] + lambda)
      }
    }
  }
  # print(paste('r: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  # rcpp code
  # from = Sys.time()
  retrofit_step3_beta(W_gk, TH_k, lambda, phi_b_gk_new, c(G,K,S));
  # print(paste('rcpp: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  expect_true(all.equal(phi_b_gk, phi_b_gk_new))
})

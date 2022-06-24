test_that("step3-alpha", {
  G = 4
  K = 3
  S = 2
  G = 1550
  K = 16
  S = 1080
  distributions = list(
    # W/H/Th from Gamma dist
    w_gk                = runif(G*K, 0, 10),
    h_ks                = runif(K*S, 0, 10),
    th_k                = runif(K, 0, 10)
  )
  phi_a_gks             = array(0, c(G,K,S))
  phi_a_gks_new         = array(0, c(G,K,S))
  lambda = runif(1, 0, 1)
  
  w_gk = array(distributions$w_gk, c(G,K))
  h_ks = array(distributions$h_ks, c(K,S))
  th_k = array(distributions$th_k, c(K))
  
  # r code
  from = Sys.time()
  for(s in 1:S){
    for(k in 1:K){
      phi_a_gks[,k,s] = ((w_gk[,k] * th_k[k]) +lambda)* h_ks[k,s]
    }
    for(v in 1:G){
      phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
    }
  }
  print(paste('r: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  # rcpp code
  from = Sys.time()
  retrofit_decomposition_step3_alpha(distributions, lambda, c(G,K,S), phi_a_gks_new);
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
  distributions = list(
    # W/H/Th from Gamma dist
    w_gk                = runif(G*K, 0, 10),
    h_ks                = runif(K*S, 0, 10),
    th_k                = runif(K, 0, 10)
  )
  phi_b_gk              = array(0, c(G,K))
  phi_b_gk_new          = array(0, c(G,K))
  lambda = runif(1, 0, 1)
  
  w_gk = array(distributions$w_gk, c(G,K))
  h_ks = array(distributions$h_ks, c(K,S))
  th_k = array(distributions$th_k, c(K))
  
  # r code
  from = Sys.time()
  for(k in 1:K){
    for(v in 1:G){
      if((w_gk[v,k]*th_k[k] + lambda)==0){
        phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lambda=0
      } else {
        phi_b_gk[v,k]=(w_gk[v,k]*th_k[k])/(w_gk[v,k]*th_k[k] + lambda)
      }
    }
  }
  print(paste('r: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  # rcpp code
  from = Sys.time()
  retrofit_decomposition_step3_beta(distributions, lambda, c(G,K,S), phi_b_gk_new);
  print(paste('rcpp: ', paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  
  expect_true(all.equal(phi_b_gk, phi_b_gk_new))
})

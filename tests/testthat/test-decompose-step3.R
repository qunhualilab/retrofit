test_that("step3-alpha-works-equally-with-Rcpp-and-R", {
  r_decompose_step3_alpha <- function(distributions, lambda, dim, phi_a_gks){
    G=dim[1]
    K=dim[2]
    S=dim[3]
    w_gk = array(distributions$w_gk, c(G,K))
    h_ks = array(distributions$h_ks, c(K,S))
    th_k = array(distributions$th_k, c(K))
    for(s in seq_len(S)){
      for(k in seq_len(K)){
        phi_a_gks[,k,s] = ((w_gk[,k] * th_k[k]) +lambda)* h_ks[k,s]
      }
      for(v in seq_len(G)){
        phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
      }
    }
    return(phi_a_gks)
  }

  G = 4
  K = 3
  S = 2
  distributions = list(
    # W/H/Th from Gamma dist
    w_gk                = runif(G*K, 0, 10),
    h_ks                = runif(K*S, 0, 10),
    th_k                = runif(K, 0, 10)
  )
  r__phi_a_gks          = array(0, c(G,K,S))
  cpp__phi_a_gks        = array(0, c(G,K,S))
  lambda = runif(1, 0, 1)
  
  # r code
  r__phi_a_gks = r_decompose_step3_alpha(distributions, lambda, c(G,K,S), r__phi_a_gks)
  
  # rcpp code
  decompose_step3_alpha(distributions, lambda, c(G,K,S), cpp__phi_a_gks)
  
  expect_true(all.equal(r__phi_a_gks, cpp__phi_a_gks))
})

test_that("step3-beta-works-equally-with-Rcpp-and-R", {
  r_decompose_step3_beta <- function(distributions, lambda, dim, phi_b_gk){
    G=dim[1]
    K=dim[2]
    S=dim[3]
    w_gk = array(distributions$w_gk, c(G,K))
    h_ks = array(distributions$h_ks, c(K,S))
    th_k = array(distributions$th_k, c(K))
    for(k in seq_len(K)){
      for(v in seq_len(G)){
        if((w_gk[v,k]*th_k[k] + lambda)==0){
          phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lambda=0
        } else {
          phi_b_gk[v,k]=(w_gk[v,k]*th_k[k])/(w_gk[v,k]*th_k[k] + lambda)
        }
      }
    }
    return(phi_b_gk)
  }
  
  G = 4
  K = 3
  S = 2
  distributions = list(
    # W/H/Th from Gamma dist
    w_gk                = runif(G*K, 0, 10),
    h_ks                = runif(K*S, 0, 10),
    th_k                = runif(K, 0, 10)
  )
  phi_b_gk              = array(0, c(G,K))
  phi_b_gk_new          = array(0, c(G,K))
  lambda = runif(1, 0, 1)
  
  # r code
  phi_b_gk = r_decompose_step3_beta(distributions, lambda, c(G,K,S), phi_b_gk)
  
  # rcpp code
  decompose_step3_beta(distributions, lambda, c(G,K,S), phi_b_gk_new);
  
  expect_true(all.equal(phi_b_gk, phi_b_gk_new))
})

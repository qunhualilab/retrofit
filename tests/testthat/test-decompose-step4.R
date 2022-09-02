test_that("step4-alpha-works-equally-with-Rcpp-and-R", {
  r_decompose_step4_alpha <- function(x, 
                                      probabilities,
                                      alpha_w_0,
                                      alpha_h_0,
                                      alpha_th_0, 
                                      dim){
    G=dim[1]
    K=dim[2]
    S=dim[3]
    alpha_w_gk = matrix(rep(0, G*K) , nrow=G, ncol=K)
    alpha_h_ks = matrix(rep(0, K*S) , nrow=K, ncol=S)
    alpha_th_k = array(rep(0, K) , c(K))
    phi_a_gks <- probabilities$phi_a_gks
    phi_b_gk <- probabilities$phi_b_gk
    for(k in 1:K){
      alpha_w_gk[,k] = alpha_w_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k]
      alpha_h_ks[k,] = alpha_h_0 + colSums(x*phi_a_gks[,k,])
      alpha_th_k[k]  = alpha_th_0+ sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
    }
    return(list(w=alpha_w_gk, h=alpha_h_ks, t=alpha_th_k))
  }
  
  G = 4
  K = 3
  S = 2
  x=matrix(runif(G*S, 0, 10),nrow=G, ncol=S)
  probabilities = list(
    phi_a_gks=array(runif(G*K*S, 0, 10), c(G,K,S)),
    phi_b_gk=array(runif(G*K, 0, 10), c(G,K))
  )
  alpha_w_0 = runif(1, 0, 1)
  alpha_h_0 = runif(1, 0, 1)
  alpha_th_0 = runif(1, 0, 1)
  
  # r code
  ret = decompose_step4_alpha(x, 
                              probabilities,
                              alpha_w_0,
                              alpha_h_0,
                              alpha_th_0, 
                              c(G,K,S))
  r__alpha_w_gk = matrix(ret$w, nrow=G, ncol=K)
  r__alpha_h_ks = matrix(ret$h, nrow=K, ncol=S)
  r__alpha_th_k = array(ret$t , c(K))
  
  # rcpp code
  ret = decompose_step4_alpha(x, 
                              probabilities,
                              alpha_w_0,
                              alpha_h_0,
                              alpha_th_0, 
                              c(G,K,S))
  cpp__alpha_w_gk = matrix(ret$w, nrow=G, ncol=K)
  cpp__alpha_h_ks = matrix(ret$h, nrow=K, ncol=S)
  cpp__alpha_th_k = array(ret$t , c(K))
  
  expect_true(all.equal(r__alpha_w_gk, cpp__alpha_w_gk))
  expect_true(all.equal(r__alpha_h_ks, cpp__alpha_h_ks))
  expect_true(all.equal(r__alpha_th_k, cpp__alpha_th_k))
})

test_that("step4-beta-works-equally-with-Rcpp-and-R", {
  r_decompose_step4_beta <- function(distributions,
                                     beta_w_0,
                                     beta_h_0,
                                     beta_th_0,
                                     lambda, 
                                     dim){
    G=dim[1]
    K=dim[2]
    S=dim[3]
    beta_w_gk=matrix(rep(0, G*K), nrow=G, ncol=K)
    beta_h_ks=matrix(rep(0, K*S), nrow=K, ncol=S)
    beta_th_k=array(rep(0, K) , c(K))
    for(k in 1:K){
      beta_w_gk[,k]= beta_w_0 + sum(h_ks[k,]*th_k[k])
      beta_h_ks[k,]= beta_h_0 + sum(w_gk[,k]*th_k[k] + lambda)
      beta_th_k[k]= beta_th_0 + sum(as.matrix(w_gk[,k]) %*% t(as.matrix(h_ks[k,])))
    }
    return(list(w=beta_w_gk, h=beta_h_ks, t=beta_th_k))
  }
  
  G = 4
  K = 3
  S = 2
  beta_w_0 = runif(1, 0, 1)
  beta_h_0 = runif(1, 0, 1)
  beta_th_0 = runif(1, 0, 1)
  w_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  h_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  th_k=runif(K, 0, 10)
  lambda = runif(1, 0, 1)
  distributions = list(
    w_gk = w_gk,
    h_ks = h_ks,
    th_k = th_k
  )
  
  # r code
  ret = r_decompose_step4_beta(distributions,
                               beta_w_0,
                               beta_h_0,
                               beta_th_0,
                               lambda, 
                               c(G,K,S))
  r__beta_w_gk = matrix(ret$w, nrow=G, ncol=K)
  r__beta_h_ks = matrix(ret$h, nrow=K, ncol=S)
  r__beta_th_k = array(ret$t , c(K))

  # rcpp code
  ret = decompose_step4_beta(distributions,
                             beta_w_0,
                             beta_h_0,
                             beta_th_0,
                             lambda, 
                             c(G,K,S))
  cpp__beta_w_gk = matrix(ret$w, nrow=G, ncol=K)
  cpp__beta_h_ks = matrix(ret$h, nrow=K, ncol=S)
  cpp__beta_th_k = array(ret$t , c(K))

  expect_true(all.equal(r__beta_w_gk, cpp__beta_w_gk))
  expect_true(all.equal(r__beta_h_ks, cpp__beta_h_ks))
  expect_true(all.equal(r__beta_th_k, cpp__beta_th_k))
})

test_that("step5-alpha-works-equally-with-Rcpp-and-R", {
  r_decompose_step5_alpha <- function(w_gk, w_gk_updated, 
                                      h_ks, h_ks_updated,
                                      th_k, th_k_updated, 
                                      rho, 
                                      dim){
    G=dim[1]
    K=dim[2]
    S=dim[3]
    alpha_w_gk = array(rep(0,G*K), c(G,K))
    alpha_h_ks = array(rep(0,K*S), c(K,S))
    alpha_th_k = array(rep(0,K), c(K))
    for(k in seq_len(K)){
      alpha_w_gk[,k]= (1-rho)*w_gk[,k] + rho*w_gk_updated[,k]
      alpha_h_ks[k,]= (1-rho)*h_ks[k,] + rho*h_ks_updated[k,]
      alpha_th_k[k]= (1-rho)*th_k[k] + rho*th_k_updated[k]
    }
    return(list(w=alpha_w_gk, h=alpha_h_ks, t=alpha_th_k))
  }
  
  G = 4
  K = 3
  S = 2
  w_gk=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  w_gk_updated=matrix(runif(G*K, 0, 10),nrow=G, ncol=K)
  h_ks=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  h_ks_updated=matrix(runif(K*S, 0, 10),nrow=K, ncol=S)
  th_k=array(runif(K, 0, 10))
  th_k_updated=array(runif(K, 0, 10))
  rho = runif(1, 0, 1)
  
  # r code
  ret = r_decompose_step5_alpha(w_gk, w_gk_updated, 
                                h_ks, h_ks_updated,
                                th_k, th_k_updated, 
                                rho, 
                                c(G,K,S))
  r__alpha_w_gk = matrix(ret$w,nrow=G, ncol=K)
  r__alpha_h_ks = matrix(ret$h,nrow=K, ncol=S)
  r__alpha_th_k = array(ret$t)
  
  # rcpp code
  cpp__alpha_w_gk <- w_gk
  cpp__alpha_h_ks <- h_ks
  cpp__alpha_th_k <- th_k
  decompose_step5(cpp__alpha_w_gk, w_gk_updated, rho)
  decompose_step5(cpp__alpha_h_ks, h_ks_updated, rho)
  decompose_step5(cpp__alpha_th_k, th_k_updated, rho)
  
  expect_true(all.equal(r__alpha_w_gk, cpp__alpha_w_gk))
  expect_true(all.equal(r__alpha_h_ks, cpp__alpha_h_ks))
  expect_true(all.equal(r__alpha_th_k, cpp__alpha_th_k))
})

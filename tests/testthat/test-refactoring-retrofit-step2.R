test_that("step2-set.seed-works-equivalently-with-Rcpp-and-R", {
  K = 4
  S = 3
  shapes=matrix(runif(K*S, 0, 1),nrow=K, ncol=S)
  rates=matrix(runif(K*S, 0, 1),nrow=K, ncol=S)
  
  set.seed(1)
  val = array(rep(0, K*S), c(K, S))
  for(s in 1:S){
    for(k in 1:K){
      val[k, s] = rgamma(1, shapes[k, s], rates[k, s])
    }
  }
  
  set.seed(1)
  val_new = array(rep(0, K*S), c(K, S))
  retrofit_decomposition_step2(shapes, rates, val_new)
  
  expect_true(all.equal(val, val_new))
})

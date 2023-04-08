test_that("step2-works-equally-with-Rcpp-and-R", {
  r_decompose_step2 <- function(shapes, rates, val){
    K <- nrow(shapes)
    S <- ncol(shapes)
    for(s in seq_len(S)){
      for(k in seq_len(K)){
        val[k, s] <- rgamma(1, shapes[k, s], rates[k, s])
      }
    }
    return(val)
  }
  
  K <- 4
  S <- 3
  shapes<-matrix(runif(K*S, 0, 1),nrow=K, ncol=S)
  rates<-matrix(runif(K*S, 0, 1),nrow=K, ncol=S)
  
  # run r function
  set.seed(1)
  val <- array(rep(0, K*S), c(K, S))
  val <- r_decompose_step2(shapes, rates, val)
  
  # run cpp function
  set.seed(1)
  val_new <- array(rep(0, K*S), c(K, S))
  decompose_step2(shapes, rates, val_new)
  
  expect_true(all.equal(val, val_new))
})

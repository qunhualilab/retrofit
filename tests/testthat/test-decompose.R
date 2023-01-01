test_that("decompose-works", {
  iterations = 10
  L = 16
  seed = 1
  
  data("TestDecomposeData")
  x   = TestDecomposeData$extra5_x
  res = retrofit::decompose(x, 
                            L=L, 
                            iterations=iterations, 
                            seed=seed,
                            verbose=TRUE)

  testthat::expect_true(all.equal(as.matrix(TestDecomposeData$results$h),  res$h,  check.attributes = FALSE))
  testthat::expect_true(all.equal(as.matrix(TestDecomposeData$results$w),  res$w,  check.attributes = FALSE))
  testthat::expect_true(all.equal(as.matrix(TestDecomposeData$results$th), res$th, check.attributes = FALSE))
})

test_that("decompose-various-x-are-accepted", {
  # matrix
  G   = 6
  S   = 8
  x   = matrix(runif(G*S, 0, 1), nrow=G, ncol=S)
  res = retrofit::decompose(x, L=16, iterations=5)
  
  # matrix all zeros
  G   = 6
  S   = 8
  x   = matrix(0,nrow=G, ncol=S)
  res = retrofit::decompose(x, L=16, iterations=5)
  
  # matrix larger dimensions
  G   = 1000
  S   = 300
  x   = matrix(runif(G*S, 0, 1), nrow=G, ncol=S)
  res = retrofit::decompose(x, L=3, iterations=5)
  
  # array 
  G   = 6
  S   = 8
  x   = array(runif(G*S, 0, 1), c(G, S))
  res = retrofit::decompose(x, L=16, iterations=5)
  
  testthat::expect_true(TRUE)
})

test_that("decompose-parameters-are-validated", {
  run_decompose <- function(args,
                            no_issue_check=FALSE,
                            warning_check=FALSE,
                            error_check=FALSE) {
    out <- tryCatch(
      {
        # don't use do.call!! it has different expression evaluation.
        decompose(x           = args$x,
                  L           = args$L,
                  iterations  = args$iterations,
                  lambda      = args$lambda,
                  seed        = args$seed,
                  alpha_w_0   = args$alpha_w_0, 
                  beta_w_0    = args$beta_w_0, 
                  alpha_h_0   = args$alpha_h_0,
                  beta_h_0    = args$beta_h_0,
                  alpha_th_0  = args$alpha_th_0,
                  beta_th_0   = args$beta_th_0,
                  kappa       = args$kappa,
                  verbose     = args$verbose)
        
        if(no_issue_check) return(TRUE)
      },
      warning=function(cond) {
        print(cond$message)
        if(warning_check) return(TRUE)
      },
      error=function(cond) {
        print(cond$message)
        if(error_check) return(TRUE)
      }
    )
    return(out)
  }
  
  base_args = list(
    x           = array(runif(4*3, 0, 1), c(4, 3)),
    L           = 16,
    iterations  = 10,
    lambda      = 0.01,
    seed        = 1,
    alpha_w_0   = 0.05, 
    beta_w_0    = 0.0001, 
    alpha_h_0   = 0.2,
    beta_h_0    = 0.2,
    alpha_th_0  = 1.25,
    beta_th_0   = 10,
    kappa       = 0.5,
    verbose     = FALSE
  )
  
  # default run
  args = base_args
  testthat::expect_true(run_decompose(args, no_issue_check = TRUE))
  
  # verbose true
  args = utils::modifyList(base_args, list(verbose=TRUE))
  testthat::expect_true(run_decompose(args, no_issue_check = TRUE))
  
  # seed can be null
  args = utils::modifyList(base_args, list(seed=NULL))
  testthat::expect_true(run_decompose(args, no_issue_check = TRUE))
  
  # warning: dimensions are large
  args = utils::modifyList(base_args, list(x = array(0, c(10000, 1000))))
  testthat::expect_true(run_decompose(args, warning_check = TRUE))
  
  # error: nulls
  # x is null
  args = utils::modifyList(base_args, list(x = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # L is null
  args = utils::modifyList(base_args, list(L = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # iterations is null
  args = utils::modifyList(base_args, list(iterations = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # lambda is null
  args = utils::modifyList(base_args, list(lambda = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # kappa is null
  args = utils::modifyList(base_args, list(kappa = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # params are null
  args = utils::modifyList(base_args, list(alpha_w_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  args = utils::modifyList(base_args, list(beta_w_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  args = utils::modifyList(base_args, list(alpha_h_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  args = utils::modifyList(base_args, list(beta_h_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  args = utils::modifyList(base_args, list(alpha_th_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  args = utils::modifyList(base_args, list(beta_th_0 = NULL))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # L is 0
  args = utils::modifyList(base_args, list(L = 0))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # iterations is 0
  args = utils::modifyList(base_args, list(iterations = 0))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  
  # kappa is 0
  args = utils::modifyList(base_args, list(kappa = 0))
  testthat::expect_true(run_decompose(args, error_check = TRUE))
})

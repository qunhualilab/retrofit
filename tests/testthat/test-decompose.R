test_that("decompose-works", {
  utils::data("testSimulationData")
  d <- testSimulationData
  
  iterations <- 10
  L <- 16
  x   <- d$extra5_x
  set.seed(1)
  res <- retrofit::decompose(x, 
                            L=L, 
                            iterations=iterations, 
                            verbose=TRUE)

  testthat::expect_true(all.equal(as.matrix(d$decompose$h),  res$h,  check.attributes = FALSE))
  testthat::expect_true(all.equal(as.matrix(d$decompose$w),  res$w,  check.attributes = FALSE))
  testthat::expect_true(all.equal(as.matrix(d$decompose$th), res$th, check.attributes = FALSE))
})

test_that("decompose-accepts-various-x", {
  # matrix
  G   <- 6
  S   <- 8
  x   <- matrix(runif(G*S, 0, 1), nrow=G, ncol=S)
  res <- retrofit::decompose(x, L=16, iterations=5)
  
  # matrix all zeros
  G   <- 6
  S   <- 8
  x   <- matrix(0,nrow=G, ncol=S)
  res <- retrofit::decompose(x, L=16, iterations=5)
  
  # matrix larger dimensions
  G   <- 1000
  S   <- 300
  x   <- matrix(runif(G*S, 0, 1), nrow=G, ncol=S)
  res <- retrofit::decompose(x, L=3, iterations=5)
  
  # array 
  G   <- 6
  S   <- 8
  x   <- array(runif(G*S, 0, 1), c(G, S))
  res <- retrofit::decompose(x, L=16, iterations=5)
  
  testthat::expect_true(TRUE)
})


test_that("decompose-validates-parameters", {
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
                  init_param  = args$init_param, 
                  lambda      = args$lambda,
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

  get_args <- function() {
    return(list(
    x           = array(runif(4*3, 0, 1), c(4, 3)),
    L           = 16,
    iterations  = 10,
    lambda      = 0.01,
    kappa       = 0.5,
    verbose     = FALSE,
    init_param  = list(
      alpha_w_0  = 0.1,
      beta_w_0   = 0.1,
      alpha_h_0  = 0.1,
      beta_h_0   = 0.1,
      alpha_th_0 = 0.1,
      beta_th_0  = 0.1,
      alpha_th_k = array(rep(0.1, 16), c(16)),
      beta_th_k  = array(rep(0.1, 16), c(16)),
      alpha_w_gk = array(rep(0.1, 4*16), c(4,16)),
      beta_w_gk  = array(rep(0.1, 4*16), c(4,16)),
      alpha_h_ks = array(rep(0.1, 16*3), c(16,3)),
      beta_h_ks  = array(rep(0.1, 16*3), c(16,3))
    )))
  }

  # seed for reproducibility
  set.seed(1)

  # default run
  args = get_args()
  testthat::expect_true(run_decompose(args, no_issue_check = TRUE))

  # verbose true
  args=get_args()
  args$verbose=TRUE
  testthat::expect_true(run_decompose(args, no_issue_check = TRUE))

  # error: nulls
  # x is null
  args=get_args()
  args$x=NULL
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # L is null
  args=get_args()
  args$L=NULL
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # iterations is null
  args=get_args()
  args$iterations=NULL
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # lambda is null
  args=get_args()
  args$lambda=NULL
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # kappa is null
  args=get_args()
  args$kappa=NULL
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # L is 0
  args=get_args()
  args$L=0
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # iterations is 0
  args=get_args()
  args$iterations=0
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # kappa is 0
  args=get_args()
  args$kappa=0
  testthat::expect_true(run_decompose(args, error_check = TRUE))

  # init_param errors
  args=get_args()
  args$init_param=0
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  args=get_args()
  args$init_param=list()
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  args=get_args()
  args$init_param=list(alpha_w_0  = "no numeric")
  testthat::expect_true(run_decompose(args, error_check = TRUE))
  args=get_args()
  args$init_param=list(alpha_w_0  = 0.1,
                       beta_w_0   = 0.1,
                       alpha_h_0  = 0.1,
                       beta_h_0   = 0.1,
                       alpha_th_0 = 0.1,
                       beta_th_0  = 0.1)
  testthat::expect_true(run_decompose(args, error_check = TRUE))
})

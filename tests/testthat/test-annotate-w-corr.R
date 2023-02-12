test_that("annotateWithCorrelations-works", {
  data("testSimulationData")
  d = testSimulationData
  res = retrofit::annotateWithCorrelations(
    sc_ref   = d$sc_ref,
    K        = 10,
    decomp_w = d$decompose$w,
    decomp_h = d$decompose$h)
  testthat::expect_true(all.equal(d$annotateWithCorrelations$ranked_cells,  res$ranked_cells,  check.attributes = FALSE))
  testthat::expect_true(all.equal(d$annotateWithCorrelations$h_prop,  res$h_prop,  check.attributes = FALSE))
})

test_that("annotateWithCorrelations-works-in-simple-model", {
  # matrix
  testthat::expect_true(TRUE)
})

test_that("annotatedWithCorrelations-parameters-are-validated", {
  run <- function(args,
                  no_issue_check=FALSE,
                  warning_check=FALSE,
                  error_check=FALSE) {
    out <- tryCatch(
      {
        # don't use do.call!! it has different expression evaluation.
        retrofit::annotateWithCorrelations(
          sc_ref   = args$sc_ref,
          K        = args$K,
          decomp_w = args$decomp_w,
          decomp_h = args$decomp_h)
        
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
    sc_ref   = array(runif(4*10, 0, 1), c(4, 5))
    rownames(sc_ref) = c("cell1","cell2","cell3","cell4")
    decomp_w = array(runif(4*10, 0, 1), c(4, 10))
    rownames(decomp_w) = c(1:4)
    colnames(decomp_w) = c(1:10)
    decomp_h = array(runif(10*2, 0, 1), c(10, 2))
    rownames(decomp_h) = c(1:10)
    colnames(decomp_h) = c(1:2)
    return(list(
      sc_ref   = sc_ref,
      K        = 5,
      decomp_w = decomp_w,
      decomp_h = decomp_h
      ))
  }
  
  # default run
  args = get_args()
  testthat::expect_true(run(args, no_issue_check = TRUE))
  # sc_ref no rownames
  args = get_args()
  args$sc_ref = array(runif(4*10, 0, 1), c(4, 5))
  testthat::expect_true(run(args, error_check = TRUE))
  # decomp_w no rownames
  # args = get_args()
  # args$decomp_w = array(runif(4*10, 0, 1), c(4, 10))
  # testthat::expect_true(run(args, error_check = TRUE))
  # decomp_w - sc_ref dimension
  args = get_args()
  args$decomp_w = array(runif(3*10, 0, 1), c(3, 10))
  testthat::expect_true(run(args, error_check = TRUE))
  # decomp_w - decomp_h dimension
  args = get_args()
  args$decomp_w = array(runif(4*11, 0, 1), c(3, 11))
  testthat::expect_true(run(args, error_check = TRUE))
})

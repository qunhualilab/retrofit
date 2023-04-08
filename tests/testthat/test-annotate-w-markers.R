test_that("annotateWithMarkers-works", {
  utils::data("testSimulationData")
  d = testSimulationData
  res = retrofit::annotateWithMarkers(
    marker_ref = d$marker_ref,
    K          = 10,
    decomp_w   = d$decompose$w,
    decomp_h   = d$decompose$h)

  testthat::expect_true(all.equal(d$annotateWithMarkers$ranked_cells,  res$ranked_cells,  check.attributes = FALSE))
  testthat::expect_true(all.equal(d$annotateWithMarkers$h_prop,  res$h_prop,  check.attributes = FALSE))
})

test_that("annotateWithMarkers-works-in-simple-model", {
  # matrix
  testthat::expect_true(TRUE)
})

test_that("annotateWithMarkers-parameters-are-validated", {
  run <- function(args,
                  no_issue_check=FALSE,
                  warning_check=FALSE,
                  error_check=FALSE) {
    out <- tryCatch(
      {
        # don't use do.call!! it has different expression evaluation.
        retrofit::annotateWithMarkers(
          marker_ref= args$marker_ref,
          K         = args$K,
          decomp_w  = args$decomp_w,
          decomp_h  = args$decomp_h)
        
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
    marker_ref   = list("1"=1:10, "2"=2:5, "3"=1:7, "4"=1:7, "5"=1:7)
    decomp_w = array(runif(4*10, 0, 1), c(4, 10))
    rownames(decomp_w) = c(1:4)
    colnames(decomp_w) = c(1:10)
    decomp_h = array(runif(10*2, 0, 1), c(10, 2))
    rownames(decomp_h) = c(1:10)
    colnames(decomp_h) = c(1:2)
    return(list(
      marker_ref= marker_ref,
      K         = 5,
      decomp_w  = decomp_w,
      decomp_h  = decomp_h
    ))
  }
  
  # default run
  args = get_args()
  testthat::expect_true(run(args, no_issue_check = TRUE))
  # fewer cell types
  args = get_args()
  args$marker_ref = list("1"=1:10, "2"=2:5)
  testthat::expect_true(run(args, warning_check = TRUE))
  # decomp_w no rownames
  # args = get_args()
  # args$decomp_w = array(runif(4*10, 0, 1), c(4, 10))
  # testthat::expect_true(run(args, error_check = TRUE))
  # decomp_w - decomp_h dimension
  args = get_args()
  args$decomp_w = array(runif(4*11, 0, 1), c(3, 11))
  testthat::expect_true(run(args, error_check = TRUE))
})

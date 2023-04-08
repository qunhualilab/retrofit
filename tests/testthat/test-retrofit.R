test_that("retrofit-works-correlations", {
  utils::data("testSimulationData")
  d = testSimulationData
  iterations = 10
  L = 16
  K = 10
  set.seed(1)
  
  res = retrofit::retrofit(
    x         = d$extra5_x,
    sc_ref    = d$sc_ref,
    iterations= iterations,
    L         = L,
    K         = K)
  
  testthat::expect_true(all.equal(d$annotateWithCorrelations$ranked_cells,  res$annotateWithCorrelations$ranked_cells,  check.attributes = FALSE))
  testthat::expect_true(all.equal(d$annotateWithCorrelations$h_prop,  res$annotateWithCorrelations$h_prop,  check.attributes = FALSE))
})

test_that("retrofit-works-markers", {
  utils::data("testSimulationData")
  d = testSimulationData
  iterations = 10
  L = 16
  K = 10
  set.seed(1)
  
  res = retrofit::retrofit(
    x         = d$extra5_x,
    marker_ref= d$marker_ref,
    iterations= iterations,
    L         = L,
    K         = K)
  
  testthat::expect_true(all.equal(d$annotateWithMarkers$ranked_cells,  res$annotateWithMarkers$ranked_cells,  check.attributes = FALSE))
  testthat::expect_true(all.equal(d$annotateWithMarkers$h_prop,  res$annotateWithMarkers$h_prop,  check.attributes = FALSE))
})


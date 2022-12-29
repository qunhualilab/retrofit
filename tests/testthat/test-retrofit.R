test_that("retrofit-works", {
  iterations = 10
  L = 16
  K = 8
  seed = 1
  
  data("TestRetrofitData")
  x           = data$colon_x
  sc_ref      = data$sc_ref
  marker_ref  = data$marker_ref

  res = retrofit::retrofit(x, 
                           sc_ref=sc_ref, 
                           marker_ref=marker_ref, 
                           iterations=iterations, 
                           L=L, 
                           K=K,
                           seed=seed,
                           verbose=TRUE)

  expect_true(all.equal(as.matrix(data$results$decompose$h),                res$decompose$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$decompose$w),                res$decompose$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$decompose$th),               res$decompose$th, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$annotateWithCorrelations$h), res$annotateWithCorrelations$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$annotateWithCorrelations$w), res$annotateWithCorrelations$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$annotateWithMarkers$h),      res$annotateWithMarkers$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(data$results$annotateWithMarkers$w),      res$annotateWithMarkers$w, check.attributes = FALSE))
})

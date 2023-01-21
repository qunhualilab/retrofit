test_that("retrofit-works", {
  iterations = 10
  L = 16
  K = 8
  
  data("testRetrofitData")
  x           = testRetrofitData$colon_x
  sc_ref      = testRetrofitData$sc_ref
  marker_ref  = testRetrofitData$marker_ref

  set.seed(1)
  res = retrofit::retrofit(x, 
                           sc_ref=sc_ref, 
                           marker_ref=marker_ref, 
                           iterations=iterations, 
                           L=L, 
                           K=K,
                           verbose=TRUE)

  expect_true(all.equal(as.matrix(testRetrofitData$results$decompose$h),                res$decompose$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$decompose$w),                res$decompose$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$decompose$th),               res$decompose$th, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$annotateWithCorrelations$h), res$annotateWithCorrelations$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$annotateWithCorrelations$w), res$annotateWithCorrelations$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$annotateWithMarkers$h),      res$annotateWithMarkers$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitData$results$annotateWithMarkers$w),      res$annotateWithMarkers$w, check.attributes = FALSE))
})

test_that("retrofit-works", {
  iterations = 10
  L = 16
  K = 8
  
  data("testRetrofitObj")
  x           = testRetrofitObj$colon_x
  sc_ref      = testRetrofitObj$sc_ref
  marker_ref  = testRetrofitObj$marker_ref

  set.seed(1)
  res = retrofit::retrofit(x, 
                           sc_ref=sc_ref, 
                           marker_ref=marker_ref, 
                           iterations=iterations, 
                           L=L, 
                           K=K,
                           verbose=TRUE)

  expect_true(all.equal(as.matrix(testRetrofitObj$results$decompose$h),                res$decompose$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$decompose$w),                res$decompose$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$decompose$th),               res$decompose$th, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$annotateWithCorrelations$h), res$annotateWithCorrelations$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$annotateWithCorrelations$w), res$annotateWithCorrelations$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$annotateWithMarkers$h),      res$annotateWithMarkers$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(testRetrofitObj$results$annotateWithMarkers$w),      res$annotateWithMarkers$w, check.attributes = FALSE))
})

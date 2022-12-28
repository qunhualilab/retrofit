test_that("main-simulation-works", {
  iterations = 10
  L = 16
  K = 8
  seed = 1
  
  data(test_main_simulation)
  x           = test_main_simulation$n10m3_x
  sc_ref      = test_main_simulation$sc_ref
  marker_ref  = test_main_simulation$marker_ref
  results     = test_main_simulation$results

  res = Retrofit(x, 
                 sc_ref=sc_ref, 
                 marker_ref=marker_ref, 
                 iterations=iterations, 
                 L=L, 
                 K=K,
                 seed=seed,
                 verbose=TRUE)

  expect_true(all.equal(as.matrix(results$decomposed$h),  res$decomposed$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$decomposed$w),  res$decomposed$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$decomposed$th), res$decomposed$th, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$annotated_correlation$h), res$annotated_correlation$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$annotated_correlation$w), res$annotated_correlation$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$annotated_marker$h), res$annotated_marker$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$annotated_marker$w), res$annotated_marker$w, check.attributes = FALSE))
})




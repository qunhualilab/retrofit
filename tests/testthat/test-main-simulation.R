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

  result = Retrofit(x, 
                ref_cor=sc_ref, 
                ref_marker=marker_ref, 
                iterations=iterations, 
                L=L, 
                K=K,
                seed=seed)
  
  expect_true(all.equal(as.matrix(results$decomp_h), result$h, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$decomp_w), result$w, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$decomp_th), result$t, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$match_cor_h), result$h_cor_map, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$match_cor_w), result$w_cor_map, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$match_marker_h), result$h_marker_map, check.attributes = FALSE))
  expect_true(all.equal(as.matrix(results$match_marker_w), result$w_marker_map, check.attributes = FALSE))
})




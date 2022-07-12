test_that("measure-error-two-norm", {
  K = 4
  S = 3
  original = array(runif(K*S, 0, 10),   c(K, S))
  inferred = array(runif(K*S, 0, 10), c(K, S))
  
  norm = decompose_compute_error_two_norm(original, inferred)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      norm_expected = norm_expected + abs(inferred[k,s]/original[k,s] - 1)
    }
  }
  
  expect_equal(norm, norm_expected)
  
  # reverse
  norm = decompose_compute_error_two_norm(inferred, original)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      norm_expected = norm_expected + abs(original[k,s]/inferred[k,s] - 1)
    }
  }
  
  expect_equal(norm, norm_expected)
  
  # include update
  original_copy = array(original,   c(K, S))
  norm = decompose_compute_and_update_error_two_norm(original, inferred)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      norm_expected = norm_expected + abs(inferred[k,s]/original_copy[k,s] - 1)
      original_copy[k,s] = inferred[k,s]
    }
  }
  
  expect_equal(norm, norm_expected)
  expect_true(all.equal(original, inferred))
})


test_that("measure-error-mat-norm", {
  K = 4
  S = 3
  original = array(runif(K*S, 0, 10),   c(K, S))
  inferred = array(runif(K*S, 0, 10), c(K, S))
  
  norm = decompose_compute_error_mat_norm(original, inferred)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      err =abs(inferred[k,s]/original[k,s] - 1)
      if(err > norm_expected){
        norm_expected = err  
      }
    }
  }
  
  expect_equal(norm, norm_expected)
  
  # reverse
  norm = decompose_compute_error_mat_norm(inferred, original)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      err =abs(original[k,s]/inferred[k,s] - 1)
      if(err > norm_expected){
        norm_expected = err  
      }
    }
  }
  
  expect_equal(norm, norm_expected)
  
  # include update
  original_copy = array(original,   c(K, S))
  norm = decompose_compute_and_update_error_mat_norm(original, inferred)
  norm_expected = 0
  for (k in 1:K){
    for (s in 1:S){
      err =abs(inferred[k,s]/original_copy[k,s] - 1)
      if(err > norm_expected){
        norm_expected = err  
      }
      original_copy[k,s] = inferred[k,s]
    }
  }
  
  expect_equal(norm, norm_expected)
  expect_true(all.equal(original, inferred))
})
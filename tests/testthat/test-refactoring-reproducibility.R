test_that("reproducibility", {
  dir = "../../results"
  file="A3_1554"
  iterations = 2
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  
  # RetrofitDecompose
  x=read.csv(in_path)
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  result = RetrofitDecompose(x, iterations=iterations, seed=1)
  
  out_h = result$h
  out_w = result$w
  out_t = result$t
  
  # RetrofitDecomposeOriginal
  x=read.csv(in_path)
  result = RetrofitDecomposeOriginal(x, iterations)
  out_original_h = result$h
  out_original_w = result$w
  out_original_t = result$t
  out_original_t = array(out_original_t)
  
  expect_true(all.equal(out_original_w, out_w))
  expect_true(all.equal(out_original_h, out_h))
  expect_true(all.equal(out_original_t, out_t))
})

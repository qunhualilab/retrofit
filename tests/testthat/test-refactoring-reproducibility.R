test_that("reproducibility", {
  # dir = "~/Research/retrofit/retrofit/results"
  dir = "../../results"
  file="A3_1554"
  paths_original = retrofit_simulation_original(dir, file, iterations=2)
  paths = retrofit_simulation(dir, file, iterations=2)
  
  w_original=read.csv(paths_original$out_w_path)
  h_original=read.csv(paths_original$out_h_path)
  t_original=read.csv(paths_original$out_t_path)
  
  w=read.csv(paths$out_w_path)
  h=read.csv(paths$out_h_path)
  t=read.csv(paths$out_t_path)
  
  expect_true(all.equal(w_original, w))
  expect_true(all.equal(h_original, h))
  expect_true(all.equal(t_original, t))
})

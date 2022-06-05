

step1__sample_rgamma <- function(shape_mat, rate_map){
  
  sampled = mapply(function(s,r) rgamma(1, shape=s, rate=r), shape_mat, rate_mat)
  return(sampled)
}

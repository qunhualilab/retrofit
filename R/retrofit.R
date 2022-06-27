#' Split a string
#'
#' @param x Spatial Transciptomics Data. Matrix(GeneExpressions, Spots)
#'
#' @return Factorization results.
#' W: Matrix(GeneExpressions, Components)
#' H: Matrix(Components, Spots)
#' $\theta$: Vector(Components)
#' @export
#'
#' 
retrofit <- function(x, iterations=4000) {
  set.seed(1)
  
  # dimensions
  G=dim(x)[1] # Gene expressions
  S=dim(x)[2] # Spots
  cell_types=8 # Cell types ('K' in the dissertation)
  K=2*cell_types # Components ('L' in the dissertation)
  dim = c(G,K,S)
    
  ### initialization
  # initial parameters
  alpha_w_0   = 0.05 
  beta_w_0    = 0.0001 
  alpha_h_0   = 0.2 
  beta_h_0    = 0.2
  alpha_th_0  = 10/cell_types
  beta_th_0   = 10
  lambda      = 0.01
  kappa       = 0.5
  # W/H/Th from Gamma dist
  dist = list(
    w_gk = array(rep(0, len=G*K), c(G,K)),
    h_ks = array(rep(0, len=K*S), c(K,S)),
    th_k = array(rep(0, len=K), c(K))
  )
  # probability variables
  prob = list(
    phi_a_gks = array(rep(0, len=G*K*S), c(G,K,S)),
    phi_b_gk  = array(rep(0, len=G*K), c(G,K))
  )
  # parameter vectors
  param = list(
    alpha_th_k  = array(runif(K,0,1)+alpha_th_0, c(K)),
    beta_th_k   = array(runif(K,0,1)+beta_th_0, c(K)),
    alpha_w_gk  = array(runif(G*K,0,0.5)+alpha_w_0, c(G,K)),
    beta_w_gk   = array(runif(G*K,0,0.005)+beta_w_0, c(G,K)),
    alpha_h_ks  = array(runif(K*S,0,0.1)+alpha_h_0, c(K,S)),
    beta_h_ks   = array(runif(K*S,0,0.5)+beta_h_0, c(K,S))
  )
  
  ## start of algorithm
  t=0
  
  while(t<iterations){
    from = Sys.time()
    t = t+1
    
    # step (1)
    rho = (t)^(-kappa)
    
    # step (2) - Update dist
    retrofit_decomposition_step2(param$alpha_h_ks, param$beta_h_ks, dist$h_ks)
    retrofit_decomposition_step2(param$alpha_th_k, param$beta_th_k, dist$th_k)
    retrofit_decomposition_step2(param$alpha_w_gk, param$beta_w_gk, dist$w_gk)
    
    # step (3) - Update prob
    retrofit_decomposition_step3_alpha(dist, lambda, dim, prob$phi_a_gks)
    retrofit_decomposition_step3_beta(dist, lambda, dim, prob$phi_b_gk)
    
    # step (4) - Get new parameters
    alpha_asterisk = retrofit_decomposition_step4_alpha(x, prob, alpha_w_0, alpha_h_0, alpha_th_0, dim)
    beta_asterisk = retrofit_decomposition_step4_beta(dist, beta_w_0, beta_h_0, beta_th_0, lambda, dim)
    
    # step (5) - Update param
    retrofit_decomposition_step5(param$alpha_w_gk, alpha_asterisk$w, rho)
    retrofit_decomposition_step5(param$alpha_h_ks, alpha_asterisk$h, rho)
    retrofit_decomposition_step5(param$alpha_th_k, alpha_asterisk$t, rho)
    retrofit_decomposition_step5(param$beta_w_gk, beta_asterisk$w, rho)
    retrofit_decomposition_step5(param$beta_h_ks, beta_asterisk$h, rho)
    retrofit_decomposition_step5(param$beta_th_k, beta_asterisk$t, rho)
    
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  }
  
  w_hat=param$alpha_w_gk/param$beta_w_gk
  h_hat=param$alpha_h_ks/param$beta_h_ks
  th_hat=param$alpha_th_k/param$beta_th_k
  result <- list(w=w_hat, h=h_hat, t=th_hat)
  return(result)
}
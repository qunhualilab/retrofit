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
  lambda       = 0.01
  kappa       = 0.5
  # memory managed vectors
  distributions = list(
    # W/H/Th from Gamma dist
    w_gk                = rep(0, len=G*K),
    h_ks                = rep(0, len=K*S),
    th_k                = rep(0, len=K)
  )
  probabilities = list(
    # probability variables
    phi_a_gks           = rep(0, len=G*K*S),
    phi_b_gk            = rep(0, len=G*K)
  )
  parameters = list(
    # parameter vectors
    alpha_th_k          = runif(K,0,1)+alpha_th_0,
    beta_th_k           = runif(K,0,1)+beta_th_0,
    alpha_w_gk          = runif(G*K,0,0.5)+alpha_w_0,
    beta_w_gk           = runif(G*K,0,0.005)+beta_w_0,
    alpha_h_ks          = runif(K*S,0,0.1)+alpha_h_0,
    beta_h_ks           = runif(K*S,0,0.5)+beta_h_0
  )
  # parameter matrices
  # alpha_TH_k  = runif(K,0,1)+alpha_th_0
  # beta_TH_k   = runif(K,0,1)+beta_th_0
  # alpha_W_gk  = runif(G*K,0,0.5)+alpha_w_0
  # beta_W_gk   = runif(G*K,0,0.005)+beta_w_0
  # alpha_H_ks  = runif(K*S,0,0.1)+alpha_h_0
  # beta_H_ks   = runif(K*S,0,0.5)+beta_h_0
  # variational parameters
  # phi_a_gks   = rep(0, len=G*K*S)
  # phi_b_gk    = rep(0, len=G*K)
  # # W/H/Th from Gamma dist
  # W_gk        = rep(0, len=G*K)
  # H_ks        = rep(0, len=K*S)
  # TH_k        = rep(0, len=K)
  
  ## start of algorithm
  t=0
  
  while(t<iterations){
    from = Sys.time()
    t = t+1
    
    # step (1)
    rho = (t)^(-kappa)
    
    # step (2) - Update distributions
    retrofit_decomposition_step2(parameters$alpha_h_ks, parameters$beta_h_ks, distributions$h_ks)
    retrofit_decomposition_step2(parameters$alpha_th_k, parameters$beta_th_k, distributions$th_k)
    retrofit_decomposition_step2(parameters$alpha_w_gk, parameters$beta_w_gk, distributions$w_gk)
    
    # step (3) - Update probabilities
    retrofit_decomposition_step3_alpha(distributions, lambda, dim, probabilities$phi_a_gks)
    retrofit_decomposition_step3_beta(distributions, lambda, dim, probabilities$phi_b_gk)
    
    # step (4) - Get new parameter values
    alpha_asterisk = retrofit_decomposition_step4_alpha(x, probabilities, alpha_w_0, alpha_h_0, alpha_th_0, dim)
    beta_asterisk = retrofit_decomposition_step4_beta(distributions, beta_w_0, beta_h_0, beta_th_0, lambda, dim)
    
    # step (5) - Update parameters
    retrofit_decomposition_step5(parameters$alpha_w_gk, alpha_asterisk$w, rho)
    retrofit_decomposition_step5(parameters$alpha_h_ks, alpha_asterisk$h, rho)
    retrofit_decomposition_step5(parameters$alpha_th_k, alpha_asterisk$t, rho)
    retrofit_decomposition_step5(parameters$beta_w_gk, beta_asterisk$w, rho)
    retrofit_decomposition_step5(parameters$beta_h_ks, beta_asterisk$h, rho)
    retrofit_decomposition_step5(parameters$beta_th_k, beta_asterisk$t, rho)
    
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  }
  
  W_hat=array(parameters$alpha_w_gk/parameters$beta_w_gk, c(G,K))
  H_hat=array(parameters$alpha_h_ks/parameters$beta_h_ks, c(K,S))
  TH_hat=array(parameters$alpha_th_k/parameters$beta_th_k, c(K))
  result <- list(w=W_hat, h=H_hat, t=TH_hat)
  return(result)
}
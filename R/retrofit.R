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
  dim = dim(x)
  G=dim[1] # Gene expressions
  S=dim[2] # Spots
  cell_types=8 # Cell types ('K' in the dissertation)
  K=2*cell_types # Components ('L' in the dissertation)
  
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
  # parameters = list(
  #   # parameter vectors
  #   alpha_th_k          = runif(K,0,1)+alpha_th_0,
  #   beta_th_k           = runif(K,0,1)+beta_th_0,
  #   alpha_w_gk          = runif(G*K,0,0.5)+alpha_w_0,
  #   beta_w_gk           = runif(G*K,0,0.005)+beta_w_0,
  #   alpha_h_ks          = runif(K*S,0,0.1)+alpha_h_0,
  #   beta_h_ks           = runif(K*S,0,0.5)+beta_h_0,
  #   alpha_th_k_asterisk = rep(K,0),
  #   beta_th_k_asterisk  = rep(K,0),
  #   alpha_w_gk_asterisk = rep(G*K,0),
  #   beta_w_gk_asterisk  = rep(G*K,0),
  #   alpha_h_ks_asterisk = rep(K*S,0),
  #   beta_h_ks_asterisk  = rep(K*S,0)
  # )
  # parameter matrices
  alpha_TH_k  = runif(K,0,1)+alpha_th_0
  beta_TH_k   = runif(K,0,1)+beta_th_0
  alpha_W_gk  = runif(G*K,0,0.5)+alpha_w_0
  beta_W_gk   = runif(G*K,0,0.005)+beta_w_0
  alpha_H_ks  = runif(K*S,0,0.1)+alpha_h_0
  beta_H_ks   = runif(K*S,0,0.5)+beta_h_0
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
    
    # step (2)
    retrofit_decomposition_step2(alpha_H_ks, beta_H_ks, distributions$h_ks)
    retrofit_decomposition_step2(alpha_TH_k, beta_TH_k, distributions$th_k)
    retrofit_decomposition_step2(alpha_W_gk, beta_W_gk, distributions$w_gk)
    
    # step (3)
    retrofit_decomposition_step3_alpha(distributions, lambda, c(G,K,S), probabilities$phi_a_gks)
    retrofit_decomposition_step3_beta(distributions, lambda, c(G,K,S), probabilities$phi_b_gk)
    
    # step (4)
    alpha_new = retrofit_step4_alpha_calculation(x, probabilities$phi_a_gks, probabilities$phi_b_gk, alpha_w_0, alpha_h_0, alpha_th_0, c(G,K,S))
    beta_new  = retrofit_step4_beta_calculation(distributions$w_gk, distributions$h_ks, distributions$th_k, beta_w_0, beta_h_0, beta_th_0, lambda, c(G,K,S))
    
    # step (5)
    alpha_W_gk = retrofit_step5_parameter_estimation(alpha_W_gk, alpha_new$w, rho)
    alpha_H_ks = retrofit_step5_parameter_estimation(alpha_H_ks, alpha_new$h, rho)
    alpha_TH_k = retrofit_step5_parameter_estimation(alpha_TH_k, alpha_new$t, rho)
    beta_W_gk  = retrofit_step5_parameter_estimation(beta_W_gk, beta_new$w, rho)
    beta_H_ks  = retrofit_step5_parameter_estimation(beta_H_ks, beta_new$h, rho)
    beta_TH_k  = retrofit_step5_parameter_estimation(beta_TH_k, beta_new$t, rho)
    
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  }
  
  W_hat=array(alpha_W_gk/beta_W_gk, c(G,K))
  H_hat=array(alpha_H_ks/beta_H_ks, c(K,S))
  TH_hat=array(alpha_TH_k/beta_TH_k, c(K))
  result <- list(w=W_hat, h=H_hat, t=TH_hat)
  return(result)
}
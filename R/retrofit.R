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
  alpha_W_0   = 0.05 
  beta_W_0    = 0.0001 
  alpha_H_0   = 0.2 
  beta_H_0    = 0.2
  alpha_TH_0  = 10/cell_types
  beta_TH_0   = 10
  lamda       = 0.01
  kappa       = 0.5
  # parameter matrices
  alpha_TH_k  = runif(K,0,1)+alpha_TH_0
  beta_TH_k   = runif(K,0,1)+beta_TH_0
  alpha_W_gk  = runif(G*K,0,0.5)+alpha_W_0
  beta_W_gk   = runif(G*K,0,0.005)+beta_W_0
  alpha_H_ks  = runif(K*S,0,0.1)+alpha_H_0
  beta_H_ks   = runif(K*S,0,0.5)+beta_H_0
  # variational parameters
  phi_a_gks   = rep(0, len=G*K*S)
  phi_b_gk    = rep(0, len=G*K)
  # W/H/Th from Gamma dist
  W_gk        = rep(0, len=G*K)
  H_ks        = rep(0, len=K*S)
  TH_k        = rep(0, len=K)
  
  ## start of algorithm
  t=0
  
  while(t<iterations){
    from = Sys.time()
    t = t+1
    
    # step (1)
    rho = (t)^(-kappa)
    
    # step (2)
    H_ks = retrofit_step2_rgamma(alpha_H_ks, beta_H_ks)
    TH_k = retrofit_step2_rgamma(alpha_TH_k, beta_TH_k)
    W_gk = retrofit_step2_rgamma(alpha_W_gk, beta_W_gk)
    
    # step (3)
    retrofit_step3_alpha(W_gk, TH_k, H_ks, lamda, phi_a_gks, c(G,K,S))
    retrofit_step3_beta(W_gk, TH_k, lamda, phi_b_gk, c(G,K,S))
    
    # step (4)
    alpha_new = retrofit_step4_alpha_calculation(x, phi_a_gks, phi_b_gk, alpha_W_0, alpha_H_0, alpha_TH_0, c(G,K,S))
    beta_new  = retrofit_step4_beta_calculation(W_gk, H_ks, TH_k, beta_W_0, beta_H_0, beta_TH_0, lamda, c(G,K,S))
    
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
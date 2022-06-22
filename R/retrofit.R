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
  # if(!is.matrix(x)) stop("x must be a matrix with shape: GeneExpressions:Spots")
  # the first column is the row names
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  # dimensions
  cell_types=8 # Cell types ('K' in the dissertation)
  G=nrow(x) # Gene expressions
  S=ncol(x) # Spots
  K=2*cell_types # Components ('L' in the dissertation)
  
  ### initialization
  # initial parameters
  alpha_W_0=0.05 
  beta_W_0=0.0001 
  alpha_H_0=0.2 
  beta_H_0=0.2
  alpha_TH_0=10/cell_types
  beta_TH_0=10
  lamda=0.01
  # parameter matrices
  alpha_TH_k=runif(K,0,1)+alpha_TH_0
  beta_TH_k=runif(K,0,1)+beta_TH_0
  alpha_W_gk=matrix(runif(G*K,0,0.5)+alpha_W_0 , nrow=G, ncol=K)
  beta_W_gk=matrix(runif(G*K,0,0.005)+beta_W_0 , nrow=G, ncol=K)
  alpha_H_ks=matrix(runif(K*S,0,0.1)+alpha_H_0 , nrow=K, ncol=S)
  beta_H_ks=matrix(runif(K*S,0,0.5)+beta_H_0 , nrow=K, ncol=S)
  # variational parameters
  phi_a_gks=array(rep(0,G*K*S), c(G,K,S))
  phi_b_gk=matrix(0, nrow=G, ncol=K)
  # W/H/Th from Gamma dist
  W_gk=matrix(0,nrow=G, ncol=K)
  H_ks=matrix(0,nrow=K, ncol=S)
  TH_k=rep(0,K)
  
  
  kappa=0.5
  
  ## start of algorithm
  t=0
  
  while(t<iterations){
    from = Sys.time()
    t=t+1
    
    # step (1)
    rho=(t)^(-kappa)
    
    # step (2)
    H_ks=matrix(retrofit_step2_rgamma(alpha_H_ks, beta_H_ks), nrow=K, ncol=S)
    TH_k=array(retrofit_step2_rgamma(alpha_TH_k, beta_TH_k), c(K)) 
    W_gk=matrix(retrofit_step2_rgamma(alpha_W_gk, beta_W_gk), nrow=G, ncol=K)
    
    # step (3)
    retrofit_step3_alpha(W_gk, TH_k, H_ks, lamda, phi_a_gks)
    retrofit_step3_beta(W_gk, TH_k, lamda, phi_b_gk)
    
    # step (4)
    alpha_updated = retrofit_step4_alpha_calculation(x, phi_a_gks, phi_b_gk, alpha_W_0, alpha_H_0, alpha_TH_0)
    beta_updated = retrofit_step4_beta_calculation(W_gk, H_ks, TH_k, beta_W_0, beta_H_0, beta_TH_0, lamda)
    
    # step (5)
    alpha_W_gk = matrix(retrofit_step5_parameter_estimation(alpha_W_gk, alpha_updated$w, rho), nrow=G, ncol=K) 
    alpha_H_ks = matrix(retrofit_step5_parameter_estimation(alpha_H_ks, alpha_updated$h, rho), nrow=K, ncol=S)
    alpha_TH_k = array(retrofit_step5_parameter_estimation(alpha_TH_k, alpha_updated$t, rho), c(K))
    beta_W_gk = matrix(retrofit_step5_parameter_estimation(beta_W_gk, beta_updated$w, rho), nrow=G, ncol=K)
    beta_H_ks = matrix(retrofit_step5_parameter_estimation(beta_H_ks, beta_updated$h, rho), nrow=K, ncol=S)
    beta_TH_k = array(retrofit_step5_parameter_estimation(beta_TH_k, beta_updated$t, rho), c(K))
    
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  }
  
  # step (6) 
  W_hat=alpha_W_gk/beta_W_gk
  H_hat=alpha_H_ks/beta_H_ks
  TH_hat=alpha_TH_k/beta_TH_k
  result <- list(w=W_hat, h=H_hat, t=TH_hat)
  return(result)
}
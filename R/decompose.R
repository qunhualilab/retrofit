#' RETROFIT decomposition algorithm
#' 
#' @description Receiving the input with 2d spatial transcriptomics matrix, 
#'  the function returns factorized {W, H, θ}.
#'  This function fulfills Structured Stochastic Variational Inference Algorithm for RETROFIT.
#'  Since exact Bayesian inference is infeasible and considering the large number of spots and genes,
#'  variational inference was adopted to approximately estimate the parameters in performant manner.
#'
#' @param x Matrix(GeneExpressions, Spots): Spatial Transciptomics Data. 
#' @param iterations integer: The number of iterations to be executed
#'
#' @return A list of decomposed vectors that contains  
#' \itemize{
#' \item w: 2d array with GeneExpressions, Components
#' \item h: 2d array with Components, Spots
#' \item th: an array with Components
#' }
#'
#'@examples
#'x=read.csv(in_path)
#'rownames(x)=x[,1]
#'x=as.matrix(x[,-1])
#'result = retrofit_decompose(x)
#'@seealso papers reference
#'@export
RetrofitDecompose <- function(x, iterations=4000) {
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
    decompose_step2(param$alpha_h_ks, param$beta_h_ks, dist$h_ks)
    decompose_step2(param$alpha_th_k, param$beta_th_k, dist$th_k)
    decompose_step2(param$alpha_w_gk, param$beta_w_gk, dist$w_gk)
    
    # step (3) - Update prob
    decompose_step3_alpha(dist, lambda, dim, prob$phi_a_gks)
    decompose_step3_beta(dist, lambda, dim, prob$phi_b_gk)
    
    # step (4) - Get new parameters
    alpha_asterisk = decompose_step4_alpha(x, prob, alpha_w_0, alpha_h_0, alpha_th_0, dim)
    beta_asterisk = decompose_step4_beta(dist, beta_w_0, beta_h_0, beta_th_0, lambda, dim)
    
    # step (5) - Update param
    decompose_step5(param$alpha_w_gk, alpha_asterisk$w, rho)
    decompose_step5(param$alpha_h_ks, alpha_asterisk$h, rho)
    decompose_step5(param$alpha_th_k, alpha_asterisk$t, rho)
    decompose_step5(param$beta_w_gk, beta_asterisk$w, rho)
    decompose_step5(param$beta_h_ks, beta_asterisk$h, rho)
    decompose_step5(param$beta_th_k, beta_asterisk$t, rho)
    
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
  }
  
  w_hat=array(param$alpha_w_gk/param$beta_w_gk, c(G,K))
  h_hat=array(param$alpha_h_ks/param$beta_h_ks, c(K,S))
  th_hat=array(param$alpha_th_k/param$beta_th_k, c(K))
  result <- list(w=w_hat, h=h_hat, t=th_hat)
  return(result)
}
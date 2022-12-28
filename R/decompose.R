#' RETROFIT decomposition algorithm
#' 
#' @description Receiving the input with 2d spatial transcriptomics matrix, 
#'  the function returns factorized {W, H, θ}.
#'  This function fulfills Structured Stochastic Variational Inference Algorithm for RETROFIT.
#'  Since exact Bayesian inference is infeasible and considering the large number of spots and genes,
#'  variational inference was adopted to approximately estimate the parameters in performant manner.
#'
#' @param x Matrix or Array with dimension (GeneExpressions, Spots). 
#'          This is the main spatial transciptomics data.
#' @param iterations integer: The number of maximum iterations to be executed
#' @param tolerance double: tolerance factor for convergence of the algorithm 
#' @param L integer: The number of components to be decomposed
#' @param alpha_w_0 double: Variational initial parameter for vector alpha_w
#' @param beta_w_0 double:  Variational initial parameter for vector beta_w
#' @param alpha_h_0 double: Variational initial parameter for vector alpha_h
#' @param beta_h_0 double:  Variational initial parameter for vector beta_h
#' @param alpha_th_0 double:Variational initial parameter for vector alpha_th
#' @param beta_th_0 double: Variational initial parameter for vector beta_th
#' @param lambda double: Background expression profile control
#' @param kappa double: Learning rate factor
#' @param seed double: Random variable seed in case the output should be deterministic
#' @param plot boolean: Plot relative errors
#'
#' @return A list of decomposed vectors that contains
#' \itemize{
#' \item w: 2d array with GeneExpressions, Components
#' \item h: 2d array with Components, Spots
#' \item th: an array with Components
#' }
#'
#'@examples
#'iterations = 10
#'L = 16
#'x=read.csv(paste("../data", "sample_x.csv", sep="/"), row.names = 1, check.names = FALSE)
#'result = RetrofitDecompose(x,
#'                           iterations=iterations,
#'                           L=L)
#'@seealso papers reference
#'@export
decompose <- function(x, 
                      L           = 16,
                      iterations  = 4000,
                      lambda      = 0.01,
                      seed        = NULL,
                      alpha_w_0   = 0.05, 
                      beta_w_0    = 0.0001, 
                      alpha_h_0   = 0.2,
                      beta_h_0    = 0.2,
                      alpha_th_0  = 1.25,
                      beta_th_0   = 10,
                      kappa       = 0.5,
                      verbose     = FALSE) {
  if (!is.null(seed)){
    set.seed(seed)  
  }
  
  # copy and 'purify' the matrix
  x_rownames = rownames(x)
  x_colnames = colnames(x)
  x = matrix(as.numeric(unlist(x)), nrow=nrow(x), ncol=ncol(x))
  
  # dimensions
  G   = dim(x)[1] # Gene expressions
  S   = dim(x)[2] # Spots
  K   = L # alias the component number
  dim = c(G,K,S)
  if (G*K*S > 4*1e+7){
    warning(paste('The dimension ( G,K,S =>', toString(dim), ') can be too 
                  large to handle in this function.'))
  }
  
  # W/H/Th from Gamma dist
  dist = list(
    w_gk = array(rep(0, len=G*K), c(G,K)),
    h_ks = array(rep(0, len=K*S), c(K,S)),
    th_k = array(rep(0, len=K),   c(K))
  )
  # probability variables
  prob = list(
    phi_a_gks = array(rep(0, len=G*K*S),c(G,K,S)),
    phi_b_gk  = array(rep(0, len=G*K),  c(G,K))
  )
  # parameter vectors
  param = list(
    alpha_th_k = array(stats::runif(K,  0,1)    +alpha_th_0,c(K)),
    beta_th_k  = array(stats::runif(K,  0,1)    +beta_th_0, c(K)),
    alpha_w_gk = array(stats::runif(G*K,0,0.5)  +alpha_w_0, c(G,K)),
    beta_w_gk  = array(stats::runif(G*K,0,0.005)+beta_w_0,  c(G,K)),
    alpha_h_ks = array(stats::runif(K*S,0,0.1)  +alpha_h_0, c(K,S)),
    beta_h_ks  = array(stats::runif(K*S,0,0.5)  +beta_h_0,  c(K,S))
  )
  
  # variables for verbose mode
  previous_param = list()
  relative_error = list()
  for (n in names(param)){
    previous_param[[n]] = rep(0, len=length(param[[n]]))
    relative_error[[n]] = c()
  }
  durations = NULL
  
  ## start of algorithm
  t=0
  while(t<iterations){
    from = Sys.time()
    t = t+1
    
    # step (1)
    rho = (t)^(-kappa)
    
    # step (2) - Sample distributions
    decompose_step2(param$alpha_h_ks, param$beta_h_ks, dist$h_ks)
    decompose_step2(param$alpha_th_k, param$beta_th_k, dist$th_k)
    decompose_step2(param$alpha_w_gk, param$beta_w_gk, dist$w_gk)
    
    # step (3) - Update probabilities
    decompose_step3_alpha(dist, lambda, dim, prob$phi_a_gks)
    decompose_step3_beta(dist, lambda, dim, prob$phi_b_gk)
    
    # step (4) - Calculate new parameters
    alpha_asterisk = decompose_step4_alpha(x, prob, alpha_w_0, alpha_h_0, alpha_th_0, dim)
    beta_asterisk  = decompose_step4_beta(dist, beta_w_0, beta_h_0, beta_th_0, lambda, dim)
    
    # step (5) - Update parameters
    decompose_step5(param$alpha_w_gk, alpha_asterisk$w, rho)
    decompose_step5(param$alpha_h_ks, alpha_asterisk$h, rho)
    decompose_step5(param$alpha_th_k, alpha_asterisk$t, rho)
    decompose_step5(param$beta_w_gk, beta_asterisk$w, rho)
    decompose_step5(param$beta_h_ks, beta_asterisk$h, rho)
    decompose_step5(param$beta_th_k, beta_asterisk$t, rho)
    
    if(verbose){
      # record error
      for (name in names(param)){
        e = decompose_compute_error_mat_norm(previous_param[[name]], param[[name]])
        decompose_update_original(previous_param[[name]], param[[name]])
        relative_error[[name]] = c(relative_error[[name]], e)
      }
      # record performance
      dur = as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs"))
      durations = c(durations, dur)
      print(paste('iteration:', t, paste0(round(dur, 3), " Seconds")))  
    }
  }
  
  w_hat=matrix(param$alpha_w_gk/param$beta_w_gk, nrow=G, ncol=K)
  h_hat=matrix(param$alpha_h_ks/param$beta_h_ks, nrow=K, ncol=S)
  th_hat=matrix(param$alpha_th_k/param$beta_th_k, nrow=K)
  
  rownames(w_hat) = x_rownames
  colnames(h_hat) = x_colnames
  
  if(!verbose){
    result <- list(w=w_hat, h=h_hat, th=th_hat)
  } else {
    # measure performance
    print(paste('Iteration mean: ', round(mean(durations), 3), " seconds", " total: ", round(sum(durations), 3), " seconds"))
    result <- list(w=w_hat, h=h_hat, th=th_hat, durations=durations, relative_error=relative_error)
  }
  
  return(result)
}
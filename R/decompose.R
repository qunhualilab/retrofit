#' RETROFIT decomposition algorithm
#' 
#' @description Receiving the input with 2d spatial transcriptomics matrix, 
#'  the function returns factorized {W, H, Theta}.
#'  This function fulfills Structured Stochastic Variational Inference Algorithm for RETROFIT.
#'  Since exact Bayesian inference is infeasible and considering the large number of spots and genes,
#'  variational inference was adopted to approximately estimate the parameters in performant manner.
#'
#' @param x matrix or array with dimension (GeneExpressions, Spots). This is the main spatial transciptomics data.
#' @param L           integer (default:16)    The number of components to be decomposed
#' @param iterations  integer (default:4000)  The number of maximum iterations to be executed
#' @param lambda      double  (default:0.01)  Background expression profile control
#' @param init_param  list                    Vatirational initial parameters
#' @details init_param specification
#' \itemize{
#'  \item alpha_w_0   double  (default:0.05)
#'  \item beta_w_0    double  (default:0.0001) 
#'  \item alpha_h_0   double  (default:0.2)   
#'  \item beta_h_0    double  (default:0.2)   
#'  \item alpha_th_0  double  (default:1.25)
#'  \item beta_th_0   double  (default:10)
#'  \item alpha_th_k  array   (default:array with dim c(K))
#'  \item beta_th_k   array   (default:array with dim c(K)),
#'  \item alpha_w_gk  array   (default:array with dim c(G,K)),
#'  \item beta_w_gk   array   (default:array with dim c(G,K)),
#'  \item alpha_h_ks  array   (default:array with dim c(K,S)),
#'  \item beta_h_ks   array   (default:array with dim c(K,S))
#' }
#' @param kappa       double  (default:0.5)   Learning rate factor
#' @param verbose     boolean (default:FALSE)
#'
#' @return A list of decomposed vectors that contains
#' \itemize{
#' \item w:             2d array with GeneExpressions, Components
#' \item h:             2d array with Components, Spots
#' \item th:            an array with Components
#' \item durations:     (verbose) durations vector (unit: second)
#' \item relative_error:(verbose) errors with pre-defined norm vector 
#' }
#' 
#'@examples
#' data("TestDecomposeData")
#' x   = TestDecomposeData$extra5_x
#' res = retrofit::decompose(x, L=16, iterations=10, verbose=TRUE)
#' W   = res$w
#' H   = res$h
#' TH  = res$th
#'@seealso papers reference
#'@export
decompose <- function(x, 
                      L           = 16,
                      iterations  = 4000,
                      init_param  = NULL,
                      lambda      = 0.01,
                      kappa       = 0.5,
                      verbose     = FALSE) {
  stopifnot(!is.null(x))
  stopifnot(!is.null(L))
  stopifnot(!is.null(lambda))
  stopifnot(!is.null(iterations))
  # stopifnot(!is.null(alpha_w_0))
  # stopifnot(!is.null(beta_w_0))
  # stopifnot(!is.null(alpha_h_0))
  # stopifnot(!is.null(beta_h_0))
  # stopifnot(!is.null(alpha_th_0))
  # stopifnot(!is.null(beta_th_0))
  stopifnot(!is.null(kappa))
  
  stopifnot(L > 0)
  stopifnot(iterations > 0)
  stopifnot(kappa > 0)
  
  # copy and 'purify' the matrix
  x_rownames = rownames(x)
  x_colnames = colnames(x)
  x = matrix(as.numeric(unlist(x)), nrow=nrow(x), ncol=ncol(x))
  # dimensions
  G   = dim(x)[1] # Gene expressions
  S   = dim(x)[2] # Spots
  K   = L # alias the component number
  dim = c(G,K,S)
  if (G*K*S > 1e+8){
    warning(paste('The dimension ( G,K,S =>', toString(dim), ') can be too large to handle in this function.'))
  }
  if(G*K*S == 0){
    stop("Dimensions cannot be 0")
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
    alpha_w_0  = 0.05, 
    beta_w_0   = 0.0001, 
    alpha_h_0  = 0.2,
    beta_h_0   = 0.2,
    alpha_th_0 = 1.25,
    beta_th_0  = 10,
    alpha_th_k = array(stats::runif(K,  0,1)    +1.25,   c(K)),
    beta_th_k  = array(stats::runif(K,  0,1)    +10,     c(K)),
    alpha_w_gk = array(stats::runif(G*K,0,0.5)  +0.05,   c(G,K)),
    beta_w_gk  = array(stats::runif(G*K,0,0.005)+0.0001, c(G,K)),
    alpha_h_ks = array(stats::runif(K*S,0,0.1)  +0.2,    c(K,S)),
    beta_h_ks  = array(stats::runif(K*S,0,0.5)  +0.2,    c(K,S))
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
    alpha_asterisk = decompose_step4_alpha(x, prob, param$alpha_w_0, param$alpha_h_0, param$alpha_th_0, dim)
    beta_asterisk  = decompose_step4_beta(dist, param$beta_w_0, param$beta_h_0, param$beta_th_0, lambda, dim)
    
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
  
  w_hat  = matrix(param$alpha_w_gk/param$beta_w_gk, nrow=G, ncol=K)
  h_hat  = matrix(param$alpha_h_ks/param$beta_h_ks, nrow=K, ncol=S)
  th_hat = matrix(param$alpha_th_k/param$beta_th_k, nrow=K)
  
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
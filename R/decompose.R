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
#' @param K integer: The number of components to be decomposed
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
#'x=read.csv(in_path, row.names=1)
#'result = retrofit_decompose(x)
#'@seealso papers reference
#'@export
RetrofitDecompose <- function(x, 
                              K           = 16,
                              alpha_w_0   = 0.05, 
                              beta_w_0    = 0.0001, 
                              alpha_h_0   = 0.2,
                              beta_h_0    = 0.2,
                              alpha_th_0  = 1.25,
                              beta_th_0   = 10,
                              lambda      = 0.01,
                              kappa       = 0.5,
                              iterations  = 4000, 
                              tolerance   = 1e-3,
                              seed        = NULL,
                              plot_convergence = FALSE) {
  if (!is.null(seed)){
    set.seed(seed)  
  }

  # copy and 'purify' the matrix
  x = matrix(as.numeric(unlist(x)), nrow=nrow(x), ncol=ncol(x))
  
  # dimensions
  G   = dim(x)[1] # Gene expressions
  S   = dim(x)[2] # Spots
  dim = c(G,K,S)
  if (G*K*S > 4*1e+7){
    warning(paste('The dimension ( G,K,S =>', toString(dim), ') might be too large to handle depending on the device.'))
    return()
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
    phi_b_gk  = array(rep(0, len=G*K),  c(G,K)),
    last_iter_phi_a_gks = array(rep(0, len=G*K*S),c(G,K,S))
  )
  # parameter vectors
  param = list(
    alpha_th_k  = array(runif(K,  0,1)  +alpha_th_0,c(K)),
    beta_th_k   = array(runif(K,  0,1)  +beta_th_0, c(K)),
    alpha_w_gk  = array(runif(G*K,0,0.5)+alpha_w_0, c(G,K)),
    beta_w_gk   = array(runif(G*K,0,0.005)+beta_w_0,c(G,K)),
    alpha_h_ks  = array(runif(K*S,0,0.1)+alpha_h_0, c(K,S)),
    beta_h_ks   = array(runif(K*S,0,0.5)+beta_h_0,  c(K,S))
  )
  
  param_last = list()
  relative_error = list()
  for (n in names(param)){
    param_last[[n]] = rep(0, len=length(param[[n]]))
    relative_error[[n]] = data.frame(iter=integer(), value=double())
  }
  
  error_window = 100
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
    
    # record error
    for (name in names(param)){
      e = decompose_compute_error_mat_norm(param_last[[name]], param[[name]])
      decompose_update_original(param_last[[name]], param[[name]])
      relative_error[[name]][t,] <- c(t, e)
    }
    # record performance
    dur = as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs"))
    durations = cbind(durations, dur)
    print(paste('iteration:', t, paste0(round(dur, 3), " Seconds")))
  }
  
  # convergence plotting
  if(plot_convergence){
    flush.console()
    plots = list()
    for (name in names(relative_error)){
      # set target data
      df <- relative_error[[name]]
      # calculate the most frequent range
      labels = -3:3
      breaks = c(0, 10^labels)
      levels = table(cut(df$value, labels= labels, breaks=breaks))
      majority_level = names(which.max(levels))[[1]]
      upto = strtoi(majority_level)
      
      for(scale in upto:(upto-1)){
        # set ceiling
        df$value[df$value>10^(scale)] = 10^(scale)
        # plot
        p  = (ggplot2::ggplot(df, ggplot2::aes(x=iter))
              + ggplot2::geom_line(ggplot2::aes(y=value), color="darkred")
              + ggplot2::ggtitle(paste0(name, " (upto 1e", ifelse(scale >= 0, "+", "-"), abs(scale), ")"))
        )
        plots[[paste0(name, scale)]] = p
      }  
    }
    gridExtra::grid.arrange(grobs=plots, ncol=2)
  }
  
  # measure performance
  iter_mean = round(mean(durations), 3)
  iter_sd = round(sd(durations), 3)
  print(paste('Iteration mean: ', iter_mean, " Seconds", ", Iteration std: ", iter_sd, " Seconds"))
  
  w_hat=array(param$alpha_w_gk/param$beta_w_gk, c(G,K))
  h_hat=array(param$alpha_h_ks/param$beta_h_ks, c(K,S))
  th_hat=array(param$alpha_th_k/param$beta_th_k, c(K))
  result <- list(w=w_hat, h=h_hat, t=th_hat)
  return(result)
}
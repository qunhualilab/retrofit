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
#' @examples
#' 
retrofit <- function(x, iterations=4000) {
  if(!is.matrix(x)) stop("x must be a matrix with shape: GeneExpressions:Spots")
  # the first column is the row names
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  # dimensions
  L=16 # Cell types
  G=nrow(x) # Gene expressions
  S=ncol(x) # Spots
  K=L/2 # Components
  
  ### initialization
  # initial parameters
  alpha_W_0=0.05 
  beta_W_0=0.0001 
  alpha_H_0=0.2 
  beta_H_0=0.2
  alpha_TH_0=10/K
  beta_TH_0=10
  lamda=0.01
  # initial matrices
  eta_theta=runif(K,0,1)+alpha_TH_0
  gam_theta=runif(K,0,1)+beta_TH_0
  eta_w=matrix(runif(G*K,0,0.5)+alpha_W_0 , nrow=G, ncol=K)
  gam_w=matrix(runif(G*K,0,0.005)+beta_W_0 , nrow=G, ncol=K)
  eta_h=matrix(runif(K*S,0,0.1)+alpha_H_0 , nrow=K, ncol=S)
  gam_h=matrix(runif(K*S,0,0.5)+beta_H_0 , nrow=K, ncol=S)
  # phi_alpha, phi_beta
  p=array(rep(0,G*K*S), c(G,K,S))
  si=matrix(0, nrow=G, ncol=K)
  # W/H/Th from Gamma dist
  Thet=rep(0,K)
  W1=matrix(0,nrow=G, ncol=K)
  H1=matrix(0,nrow=K, ncol=S)
  
  
  kappa=0.5
  
  ## start of algorithm
  t=1
  
  while(t<iterations){ 
    # step (1)
    rho=(t)^(-kappa)
    
    if(t %% 100 == 0){
      print(paste("Iteration",t))
    }
    
    # simulating from q(.)
    # step (2)
    for(k in 1:K){
      for(s in 1:S){
        H1[k,s]=rgamma(1, shape=eta_h[k,s], rate=gam_h[k,s])
      }
      Thet[k]=rgamma(1, shape=eta_theta[k], rate=gam_theta[k])
      for(v in 1:G){
        W1[v,k]=rgamma(1, shape=eta_w[v,k], rate=gam_w[v,k])
        
        # step (3) - phi_beta
        if((W1[v,k]*Thet[k] + lamda)==0){
          si[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
        } else{
          si[v,k]=(W1[v,k]*Thet[k])/(W1[v,k]*Thet[k] + lamda)
        }
      }
    }
    
    # step (3) - phi_alpha
    for(s in 1:S){
      for(k in 1:K){
        p[,k,s]=((W1[,k] * Thet[k]) +lamda)* H1[k,s]
      }
      for(v in 1:G){
        p[v,,s]=p[v,,s]/sum(p[v,,s])
      }
    }
    
    # for clear assignment
    eta_theta0=eta_theta
    gam_theta0=gam_theta
    eta_h0=eta_h
    gam_h0=gam_h
    eta_w0=eta_w
    gam_w0=gam_w
    
    # step (4) + (5)
    for(k in 1:K){
      eta_w[,k]= (1-rho)*eta_w0[,k] + rho*(alpha_W_0 + rowSums(x*p[,k,])*si[,k])
      
      gam_w[,k]= (1-rho)*gam_w0[,k] + rho*(beta_W_0 + sum(H1[k,]*Thet[k]))
      
      eta_h[k,]= (1- rho)*eta_h0[k,] + rho*(alpha_H_0 + colSums(x*p[,k,]))
      
      gam_h[k,]= (1- rho)*gam_h0[k,] + rho*(beta_H_0 + sum(W1[,k]*Thet[k] + lamda))
      
      eta_theta[k]= (1-rho)*eta_theta0[k] + rho*(alpha_TH_0 + sum(rowSums(x*p[,k,])*si[,k]))
      
      gam_theta[k]= (1-rho)*gam_theta0[k] + rho*(beta_TH_0 + sum(as.matrix(W1[,k]) %*% 
                                                           t(as.matrix(H1[k,]))))
      
    }
    
    t=t+1 
    if(sum(is.nan(p))>0){
      break
    }
  }
  
  # step (6) 
  W_hat=eta_w/gam_w
  H_hat=eta_h/gam_h
  Theta_hat=eta_theta/gam_theta
  result <- list(w=W_hat, h=H_hat, t=Theta_hat)
  return(result)
}
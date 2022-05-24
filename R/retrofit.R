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
  # parameter matrices
  alpha_W_gk=matrix(runif(G*K,0,0.5)+alpha_W_0 , nrow=G, ncol=K)
  beta_W_gk=matrix(runif(G*K,0,0.005)+beta_W_0 , nrow=G, ncol=K)
  alpha_H_ks=matrix(runif(K*S,0,0.1)+alpha_H_0 , nrow=K, ncol=S)
  beta_H_ks=matrix(runif(K*S,0,0.5)+beta_H_0 , nrow=K, ncol=S)
  alpha_TH_k=runif(K,0,1)+alpha_TH_0
  beta_TH_k=runif(K,0,1)+beta_TH_0
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
        H1[k,s]=rgamma(1, shape=alpha_H_ks[k,s], rate=beta_H_ks[k,s])
      }
      Thet[k]=rgamma(1, shape=alpha_TH_k[k], rate=beta_TH_k[k])
      for(v in 1:G){
        W1[v,k]=rgamma(1, shape=alpha_W_gk[v,k], rate=beta_W_gk[v,k])
        
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
    alpha_TH_k0=alpha_TH_k
    beta_TH_k0=beta_TH_k
    alpha_H_ks0=alpha_H_ks
    beta_H_ks0=beta_H_ks
    alpha_W_gk0=alpha_W_gk
    beta_W_gk0=beta_W_gk
    
    # step (4) + (5)
    for(k in 1:K){
      alpha_W_gk[,k]= (1-rho)*alpha_W_gk0[,k] + rho*(alpha_W_0 + rowSums(x*p[,k,])*si[,k])
      
      beta_W_gk[,k]= (1-rho)*beta_W_gk0[,k] + rho*(beta_W_0 + sum(H1[k,]*Thet[k]))
      
      alpha_H_ks[k,]= (1- rho)*alpha_H_ks0[k,] + rho*(alpha_H_0 + colSums(x*p[,k,]))
      
      beta_H_ks[k,]= (1- rho)*beta_H_ks0[k,] + rho*(beta_H_0 + sum(W1[,k]*Thet[k] + lamda))
      
      alpha_TH_k[k]= (1-rho)*alpha_TH_k0[k] + rho*(alpha_TH_0 + sum(rowSums(x*p[,k,])*si[,k]))
      
      beta_TH_k[k]= (1-rho)*beta_TH_k0[k] + rho*(beta_TH_0 + sum(as.matrix(W1[,k]) %*% 
                                                           t(as.matrix(H1[k,]))))
      
    }
    
    t=t+1 
    if(sum(is.nan(p))>0){
      break
    }
  }
  
  # step (6) 
  W_hat=alpha_W_gk/beta_W_gk
  H_hat=alpha_H_ks/beta_H_ks
  Theta_hat=alpha_TH_k/beta_TH_k
  result <- list(w=W_hat, h=H_hat, t=Theta_hat)
  return(result)
}
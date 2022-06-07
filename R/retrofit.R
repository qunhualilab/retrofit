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
  # variational parameters
  phi_a_gks=array(rep(0,G*K*S), c(G,K,S))
  phi_b_gk=matrix(0, nrow=G, ncol=K)
  # W/H/Th from Gamma dist
  W_gk=matrix(0,nrow=G, ncol=K)
  H_ks=matrix(0,nrow=K, ncol=S)
  TH_k=rep(0,K)
  
  
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
        H_ks[k,s]=rgamma(1, shape=alpha_H_ks[k,s], rate=beta_H_ks[k,s])
      }
      TH_k[k]=rgamma(1, shape=alpha_TH_k[k], rate=beta_TH_k[k])
      for(v in 1:G){
        W_gk[v,k]=rgamma(1, shape=alpha_W_gk[v,k], rate=beta_W_gk[v,k])
        
        # step (3) - phi_beta
        if((W_gk[v,k]*TH_k[k] + lamda)==0){
          phi_b_gk[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
        } else{
          phi_b_gk[v,k]=(W_gk[v,k]*TH_k[k])/(W_gk[v,k]*TH_k[k] + lamda)
        }
      }
    }
    
    # step (3) - phi_alpha
    for(s in 1:S){
      for(k in 1:K){
        phi_a_gks[,k,s]=((W_gk[,k] * TH_k[k]) +lamda)* H_ks[k,s]
      }
      for(v in 1:G){
        phi_a_gks[v,,s]=phi_a_gks[v,,s]/sum(phi_a_gks[v,,s])
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
      alpha_W_gk[,k]= (1-rho)*alpha_W_gk0[,k] + rho*(alpha_W_0 + rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k])
      
      beta_W_gk[,k]= (1-rho)*beta_W_gk0[,k] + rho*(beta_W_0 + sum(H_ks[k,]*TH_k[k]))
      
      alpha_H_ks[k,]= (1- rho)*alpha_H_ks0[k,] + rho*(alpha_H_0 + colSums(x*phi_a_gks[,k,]))
      
      beta_H_ks[k,]= (1- rho)*beta_H_ks0[k,] + rho*(beta_H_0 + sum(W_gk[,k]*TH_k[k] + lamda))
      
      alpha_TH_k[k]= (1-rho)*alpha_TH_k0[k] + rho*(alpha_TH_0 + sum(rowSums(x*phi_a_gks[,k,])*phi_b_gk[,k]))
      
      beta_TH_k[k]= (1-rho)*beta_TH_k0[k] + rho*(beta_TH_0 + sum(as.matrix(W_gk[,k]) %*% 
                                                           t(as.matrix(H_ks[k,]))))
      
    }
    
    t=t+1 
    if(sum(is.nan(phi_a_gks))>0){
      break
    }
  }
  
  # step (6) 
  W_hat=alpha_W_gk/beta_W_gk
  H_hat=alpha_H_ks/beta_H_ks
  TH_hat=alpha_TH_k/beta_TH_k
  result <- list(w=W_hat, h=H_hat, t=TH_hat)
  return(result)
}
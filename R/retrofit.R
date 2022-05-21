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
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  K=8
  G=nrow(x)
  S=ncol(x)
  
  ### initialization
  L=2*K
  e=0.05 
  f=0.0001 
  g=0.2 
  h=0.2
  m=10/K
  n=10
  lamda=0.01
  
  eta_theta=runif(L,0,1)+m
  gam_theta=runif(L,0,1)+n
  eta_w=matrix(runif(G*L,0,0.5)+e , nrow=G, ncol=L)
  gam_w=matrix(runif(G*L,0,0.005)+f , nrow=G, ncol=L)
  eta_h=matrix(runif(L*S,0,0.1)+g , nrow=L, ncol=S)
  gam_h=matrix(runif(L*S,0,0.5)+h , nrow=L, ncol=S)
  
  p=array(rep(0,G*L*S), c(G,L,S))
  si=matrix(0, nrow=G, ncol=L)
  
  Thet=rep(0,L)
  W1=matrix(0,nrow=G, ncol=L)
  H1=matrix(0,nrow=L, ncol=S)
  
  
  kappa=0.5
  
  ## start of algorithm
  t=1
  
  while(t<iterations){ 
    rho=(t)^(-kappa) #step size
    
    if(t %% 100 == 0){
      print(paste("Iteration",t))
    }
    
    # simulating from q(.)
    for(k in 1:L){
      for(s in 1:S){
        H1[k,s]=rgamma(1, shape=eta_h[k,s], rate=gam_h[k,s])
      }
      Thet[k]=rgamma(1, shape=eta_theta[k], rate=gam_theta[k])
      for(v in 1:G){
        W1[v,k]=rgamma(1, shape=eta_w[v,k], rate=gam_w[v,k])
        if((W1[v,k]*Thet[k] + lamda)==0){
          si[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
        } else{
          si[v,k]=(W1[v,k]*Thet[k])/(W1[v,k]*Thet[k] + lamda)
        }
      }
    }
    
    for(s in 1:S){
      for(k in 1:L){
        p[,k,s]=((W1[,k] * Thet[k]) +lamda)* H1[k,s]
      }
      for(v in 1:G){
        p[v,,s]=p[v,,s]/sum(p[v,,s])
      }
    }
    
    eta_theta0=eta_theta
    gam_theta0=gam_theta
    eta_h0=eta_h
    gam_h0=gam_h
    eta_w0=eta_w
    gam_w0=gam_w
    
    for(k in 1:L){
      eta_w[,k]= (1-rho)*eta_w0[,k] + rho*(e + rowSums(x*p[,k,])*si[,k])
      
      gam_w[,k]= (1-rho)*gam_w0[,k] + rho*(f + sum(H1[k,]*Thet[k]))
      
      eta_h[k,]= (1- rho)*eta_h0[k,] + rho*(g + colSums(x*p[,k,]))
      
      gam_h[k,]= (1- rho)*gam_h0[k,] + rho*(h + sum(W1[,k]*Thet[k] + lamda))
      
      eta_theta[k]= (1-rho)*eta_theta0[k] + rho*(m + sum(rowSums(x*p[,k,])*si[,k]))
      
      gam_theta[k]= (1-rho)*gam_theta0[k] + rho*(n + sum(as.matrix(W1[,k]) %*% 
                                                           t(as.matrix(H1[k,]))))
      
    }
    
    t=t+1 
    if(sum(is.nan(p))>0){
      break
    }
  }
  
  
  W_hat=eta_w/gam_w
  H_hat=eta_h/gam_h
  Theta_hat=eta_theta/gam_theta
  result <- list(w=W_hat, h=H_hat, t=Theta_hat)
  return(result)
}
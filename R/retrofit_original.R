# rm(list=ls())
# 
# invokeRestart("abort")

### Read synthetic data
retrofit_original <- function(X, iterations=4000) {
  set.seed(1)
  # X=read.csv("A3_counts_G=1554.csv") # change this to other ST data
  rownames(X)=X[,1]
  X=as.matrix(X[,-1])
  
  K=8
  G=nrow(X)
  S=ncol(X)
  
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
  
  while(t<=iterations){
    from = Sys.time()
    
    rho=(t)^(-kappa) #step size
  
    # simulating from q(.)
    # for(k in 1:L){
    #   for(s in 1:S){
    #     H1[k,s]=rgamma(1, shape=eta_h[k,s], rate=gam_h[k,s])
    #   }
    #   Thet[k]=rgamma(1, shape=eta_theta[k], rate=gam_theta[k])
    #   for(v in 1:G){
    #     W1[v,k]=rgamma(1, shape=eta_w[v,k], rate=gam_w[v,k])
    #     if((W1[v,k]*Thet[k] + lamda)==0){
    #       si[v,k]=1 ## to avoid numerical error of 0/0 when lamda=0
    #     } else{
    #       si[v,k]=(W1[v,k]*Thet[k])/(W1[v,k]*Thet[k] + lamda)
    #     }
    #   }
    # }
    
    # simulating from q(.) - modified to match the usage of set.seed().
    for(s in 1:S){
      for(k in 1:L){
        H1[k,s]=rgamma(1, shape=eta_h[k,s], rate=gam_h[k,s])
      }
    }
    
    for(k in 1:L){
      Thet[k]=rgamma(1, shape=eta_theta[k], rate=gam_theta[k])
    }
    
    for(k in 1:L){
      for(v in 1:G){
        W1[v,k]=rgamma(1, shape=eta_w[v,k], rate=gam_w[v,k])
      }
    }

    for(k in 1:L){
      for(v in 1:G){
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
      eta_w[,k]= (1-rho)*eta_w0[,k] + rho*(e + rowSums(X*p[,k,])*si[,k])

      gam_w[,k]= (1-rho)*gam_w0[,k] + rho*(f + sum(H1[k,]*Thet[k]))

      eta_h[k,]= (1- rho)*eta_h0[k,] + rho*(g + colSums(X*p[,k,]))

      gam_h[k,]= (1- rho)*gam_h0[k,] + rho*(h + sum(W1[,k]*Thet[k] + lamda))

      eta_theta[k]= (1-rho)*eta_theta0[k] + rho*(m + sum(rowSums(X*p[,k,])*si[,k]))

      gam_theta[k]= (1-rho)*gam_theta0[k] + rho*(n + sum(as.matrix(W1[,k]) %*%
                                                           t(as.matrix(H1[k,]))))
  
    }
  
    print(paste('iteration:', t, paste0(round(as.numeric(difftime(time1 = Sys.time(), time2 = from, units = "secs")), 3), " Seconds")))
    
    t=t+1
    if(sum(is.nan(p))>0){
      break
    }
  }
  
  
  W_hat=eta_w/gam_w
  H_hat=eta_h/gam_h
  Theta_hat=eta_theta/gam_theta
  # write.csv(H_hat,"A3_H_hat_G=1554.csv")
  # write.csv(W_hat,"A3_W_hat_G=1554.csv")
  # write.csv(Theta_hat,"A3_Theta_hat_G=1554.csv")
  result <- list(w=W_hat, h=H_hat, t=Theta_hat)
}

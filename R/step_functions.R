alpha_H_ks=matrix(runif(2*4,0,0.1)+0.2 , nrow=2, ncol=4)
print(alpha_H_ks)

product <- function(A, B, dim){
  l = Apply(list(A, B), target_dims=list(NULL,NULL), function(a, b) a*b)
  arr = array(unlist(l), dim=dim)
  return(arr)
}

# step2: rgamma sampling
sample_rgamma <- function(shape_mat, rate_mat){
  
  sampled = mapply(function(s,r) rgamma(1, shape=s, rate=r), shape_mat, rate_mat)
  sampled_mat = matrix(sampled)
  return(sampled)
}

lambda = 0.1
G = 3
K = 2
S = 4
phi_a_gks=array(rep(0,G*K*S), c(G,K,S))
phi_b_gk=matrix(0, nrow=G, ncol=K)
# W/H/Th from Gamma dist
W_gk=matrix(c(1,2,3,4,5,6),nrow=G, ncol=K)
H_ks=matrix(c(10,20,30,40,50,60,70,80),nrow=K, ncol=S)
TH_k=array(c(100,200))

# step3: phi computing
# alpha
W_TH_gk = product(W_gk, TH_k, c(G,K)) + lambda
W_TH_H_gks = product(W_TH_gk, H_ks, c(G,K,S))
W_TH_H_gks_sumk = apply(W_TH_H_gks, c(1,3), sum)
W_TH_H_gks_sumk_dup = aperm(array(rep(W_TH_H_gks_sumk, K), c(G,S,K)), c(1,3,2))
phi_a_gks = W_TH_H_gks/W_TH_H_gks_sumk_dup

# beta
W_TH_gk = product(W_gk, TH_k, c(G,K))
phi_b_gk = Apply(list(W_TH_gk), target_dims = list(NULL), function(x) ifelse((x + lambda)==0, 1, x/(x + lambda)))
phi_b_gk = array(unlist(phi_b_gk), dim=c(G,K))
# print(phi_b_gk)


source("../EFM algorithm.R")
source("../QFM algorithm.R")

Z = readRDS("CRSP(240months).rds")

###############################################################
########## Mean FM factor number ################################
###############################################################
Te <- nrow(Z)
N <- ncol(Z)
icp<-function(r){
  Ze<-scale(Z,scale=FALSE)
  svd = svd(Ze) 
  lambda = t(as.matrix(sqrt(N)*t(svd$v)[1:r,])) 
  Factor = (Ze %*% lambda)/N 
  IC = log(sum((Ze-Factor%*%t(lambda))^2)/(N*Te))+r*(N+Te)/(N*Te)*log(min(N,Te)) #icp2
  return(IC)
}
ic=c()
for(r in 2:10){
  ic[r]=icp(r) #r=6
}
which.min(ic)

###############################################################
########## QFM factor number ################################
###############################################################
tau.vec= c(0.001, 0.005, seq(0.01,0.1,0.01))
for(tau in tau.vec){
  set.seed(1)
  QFM0 = qfm_est(scale(Z,scale=FALSE), 15, tol=1e-5, tau=tau)
  lmat = QFM0$lmat
  dum = t(lmat) %*% lmat / N
  rhat.q <- sum(eigen(dum)$values > min(N,Te)^(-1/3)* eigen(dum)$values[1]) #log(N*Te/(N+Te))*(N+Te)/N/Te 
  print(rhat.q) # =1
}

###############################################################
########## EFM factor number ##################################
###############################################################
N = ncol(Z)
Te = nrow(Z)
tol=1e-5

expectile_loss <- function(u, tau) {
  w <- ifelse(u < 0, 1 - tau, tau)
  w * u^2
}

expectile_emp <- function(x, tau = tau, tol = 1e-8, maxit = 1000) {
  e <- mean(x)
  for (i in 1:maxit) {
    w <- ifelse(x < e, 1 - tau, tau)
    e_new <- sum(w * x) / sum(w)
    if (abs(e_new - e) < tol) break
    e <- e_new
  }
  e
}

icp <- function(r,tau){
  Ze = scale(Z, scale=FALSE)
  EFM0 = efm_est(Ze, r=r, tau=tau)
  lambda = EFM0$lmat
  Factor = EFM0$fmat
  exp=expectile_emp(Ze, tau=tau)
  penalty = (min(sqrt(N),sqrt(Te)))^(-2/3)*mean(expectile_loss(Ze-exp, tau=tau))*1/log(min(sqrt(N),sqrt(Te)))^2
  IC = mean(expectile_loss(Ze-Factor%*%t(lambda), tau=tau)) + r*penalty
  return(IC)
}

rhat.e = c() 

for(i in seq_along(tau.vec)){
  tau=tau.vec[i]
  IC = c()
  for(r in 1:8){
    set.seed(0)
    IC[r] = icp(r,tau) 
  }
  rhat.e[i] = which.min(IC)
  print(rhat.e[i])
}

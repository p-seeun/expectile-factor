source("../EFM algorithm.R")
Z = readRDS("FRED-MD.rds")  

###############################################################
########## Factor number ##################################
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

icp<-function(r,tau){
  Ze = scale(Z, scale=FALSE)
  EFM0 = efm_est(Ze, r=r, tau=tau) 
  lambda = EFM0$lmat
  Factor = EFM0$fmat
  exp=expectile_emp(Ze, tau=tau)
  penalty = (min(sqrt(N),sqrt(Te)))^(-2/3)*mean(expectile_loss(Ze-exp, tau=tau))*1/log(min(sqrt(N),sqrt(Te)))^2
  IC = mean(expectile_loss(Ze-Factor%*%t(lambda), tau=tau)) + r*penalty
  return(IC)
}

tau=0.5
IC = c()
for(r in 1:8){
  set.seed(1)
  IC[r] = icp(r,tau) 
}
which.min(IC) 

source("../EFM algorithm.R")
source("../QFM algorithm.R")

library(sn)

rep = 100
dist = "normal"
if(dist == "normal"){
  part_moment <- function(a) {
    integrand <- function(x) {
      x * dnorm(x)
    }
    integrate(integrand, lower = -Inf, upper = a)$value
  }
  h <- function(alpha){
    q_alpha = qnorm(alpha)
    G = part_moment(q_alpha)
    mean = 0
    h_alpha = (-alpha * q_alpha + G) / (2*G-mean + (1-2*alpha) * q_alpha)
    return(h_alpha)
  }
}

invert_h <- function(h0, lower = 1e-6, upper = 1 - 1e-6) {
  f <- function(alpha) {
    h(alpha) - h0
  }
  uniroot(f, lower = lower, upper = upper)$root
}

invert_h(0.999)


rsquare <- function(a, B){
  lm<- summary(lm(a ~ B))
  return(lm$adj.r.squared)
}

rmse <- function(A, B){
  return(sqrt(mean((A-B)^2)))
}

trace_mat <- function(A) sum(diag(A))

D <- function(O1, O2){
  q1 <- ncol(O1); q2 <- ncol(O2)
  q <- max(q1, q2)
  val <- 1 - (trace_mat(O1 %*% t(O1) %*% O2 %*% t(O2))/nrow(O1)/nrow(O2) / q)
  return(sqrt(val))
}

## QFMs
QFM0.95 = list()
QFM0.99 = list()
QFM0.995 = list()
QFM0.999 = list()

EFM0.95 = list()
EFM0.99 = list()
EFM0.995 = list()
EFM0.999 = list()

#n=1 or n=2 (if n=1, N=T=100 and if n=2, N=T=200)
n=1

# Set the parameters
N = 100*n
Te = 100*n
tol = 1e-5
b = c(0.8, 0.5, 0.2)

# Initialize vectors to store R2 results
R2_QFM0.95 = matrix(NA, nrow = rep, ncol = 3)
R2_QFM0.99 = matrix(NA, nrow = rep, ncol = 3)
R2_QFM0.995 = matrix(NA, nrow = rep, ncol = 3)
R2_QFM0.999 = matrix(NA, nrow = rep, ncol = 3)

R2_EFM0.95 = matrix(NA, nrow = rep, ncol = 3)
R2_EFM0.99 = matrix(NA, nrow = rep, ncol = 3)
R2_EFM0.995 = matrix(NA, nrow = rep, ncol = 3)
R2_EFM0.999 = matrix(NA, nrow = rep, ncol = 3)

# Initialize vectors to store RMSE results
recons_QFM0.95 = c()
recons_QFM0.99 = c()
recons_QFM0.995 = c()
recons_QFM0.999 = c()

recons_EFM0.95 = c()
recons_EFM0.99 = c()
recons_EFM0.995 = c()
recons_EFM0.999 = c()

# Initialize vectors to store AVE_FS results

D_QFM0.95 = c()
D_QFM0.99 = c()
D_QFM0.995 = c()
D_QFM0.999 = c()

D_EFM0.95 = c()
D_EFM0.99 = c()
D_EFM0.995 = c()
D_EFM0.999 = c()

q0.95=invert_h(0.95)
q0.99=invert_h(0.99)
q0.995=invert_h(0.995)
q0.999=invert_h(0.999)

for (s in 1:rep) {
  set.seed(s)
  
  # Generate data
  F <- matrix(0, nrow = Te, ncol = 3)
  u <- matrix(rnorm(Te * 3), nrow = Te, ncol = 3)
  F[1, ] <- u[1, ]
  for (h in 2:Te) {
    F[h, 1] <- b[1] * F[h - 1, 1] + u[h, 1]
    F[h, 2] <- b[2] * F[h - 1, 2] + u[h, 2]
    F[h, 3] <- b[3] * F[h - 1, 3] + u[h, 3]
  }
  Lambda <- matrix(rnorm(N * 3), nrow = N, ncol = 3)

  if(dist =="normal"){
    e =  matrix(rnorm(Te * N), nrow = Te, ncol = N)
  }
  
  X <- F %*% t(Lambda) + e
  
  # Estimate factors
  qfm0.95 = qfm_est(X, 4, tol, tau=q0.95)
  qfm0.99 = qfm_est(X, 4, tol, tau=q0.99)
  qfm0.995 = qfm_est(X, 4, tol, tau=q0.995)
  qfm0.999 = qfm_est(X, 4, tol, tau=q0.999)
  
  efm0.95 = efm_est(X, 4, tol, tau=0.95)
  efm0.99 = efm_est(X, 4, tol, tau=0.99)
  efm0.995 = efm_est(X, 4, tol, tau=0.995)
  efm0.999 = efm_est(X, 4, tol, tau=0.999)
  
  #Store recons results
  recons_QFM0.95[s] = rmse(F %*% t(Lambda) + qnorm(q0.95) , qfm0.95$fmat %*% t(qfm0.95$lmat))
  recons_QFM0.99[s] = rmse(F %*% t(Lambda) + qnorm(q0.99) , qfm0.99$fmat %*% t(qfm0.99$lmat))
  recons_QFM0.995[s] = rmse(F %*% t(Lambda) + qnorm(q0.995) , qfm0.995$fmat %*% t(qfm0.995$lmat))
  recons_QFM0.999[s] = rmse(F %*% t(Lambda) + qnorm(q0.999) , qfm0.999$fmat %*% t(qfm0.999$lmat))
  
  recons_EFM0.95[s] = rmse(F %*% t(Lambda) + qnorm(q0.95) , efm0.95$fmat %*% t(efm0.95$lmat))
  recons_EFM0.99[s] = rmse(F %*% t(Lambda) + qnorm(q0.99) , efm0.99$fmat %*% t(efm0.99$lmat))
  recons_EFM0.995[s] = rmse(F %*% t(Lambda) + qnorm(q0.995) , efm0.995$fmat %*% t(efm0.995$lmat))
  recons_EFM0.999[s] = rmse(F %*% t(Lambda) + qnorm(q0.999) , efm0.999$fmat %*% t(efm0.999$lmat))
  
  # Store R2 results
  R2_QFM0.95[s, ] <- apply(F, 2, function(x) rsquare(x, qfm0.95$fmat))
  R2_QFM0.99[s, ] <- apply(F, 2, function(x) rsquare(x, qfm0.99$fmat))
  R2_QFM0.995[s, ] <- apply(F, 2, function(x) rsquare(x, qfm0.995$fmat))
  R2_QFM0.999[s, ] <- apply(F, 2, function(x) rsquare(x, qfm0.999$fmat))
  
  R2_EFM0.95[s, ] <- apply(F, 2, function(x) rsquare(x, efm0.95$fmat))
  R2_EFM0.99[s, ] <- apply(F, 2, function(x) rsquare(x, efm0.99$fmat))
  R2_EFM0.995[s, ] <- apply(F, 2, function(x) rsquare(x, efm0.995$fmat))
  R2_EFM0.999[s, ] <- apply(F, 2, function(x) rsquare(x, efm0.999$fmat))
  
  #Normalize F
  F <- cbind(1,F)
  S <- crossprod(F) / Te          # (r x r) = t(F)%*%F / T
  R <- chol(S)                     # S = t(R) %*% R  (upper-triangular)
  
  Fhat <- F %*% solve(R) 
  
  D_QFM0.95[s] = D(Fhat, qfm0.95$fmat)
  D_QFM0.99[s] = D(Fhat, qfm0.99$fmat)
  D_QFM0.995[s] = D(Fhat, qfm0.995$fmat)
  D_QFM0.999[s] = D(Fhat, qfm0.999$fmat)
  
  D_EFM0.95[s] = D(Fhat, efm0.95$fmat)
  D_EFM0.99[s] = D(Fhat, efm0.99$fmat)
  D_EFM0.995[s] = D(Fhat, efm0.995$fmat)
  D_EFM0.999[s] = D(Fhat, efm0.999$fmat)
  
  print(s)
}

Mean_R2_QFM0.95  <- colMeans(R2_QFM0.95); sd_R2_QFM0.95  <- apply(R2_QFM0.95, 2, sd)
Mean_R2_QFM0.99  <- colMeans(R2_QFM0.99); sd_R2_QFM0.99  <- apply(R2_QFM0.99, 2, sd)
Mean_R2_QFM0.995 <- colMeans(R2_QFM0.995); sd_R2_QFM0.995  <- apply(R2_QFM0.995, 2, sd)
Mean_R2_QFM0.999 <- colMeans(R2_QFM0.999); sd_R2_QFM0.999  <- apply(R2_QFM0.999, 2, sd)

Mean_R2_EFM0.95  <- colMeans(R2_EFM0.95); sd_R2_EFM0.95  <- apply(R2_EFM0.95, 2, sd)
Mean_R2_EFM0.99  <- colMeans(R2_EFM0.99); sd_R2_EFM0.99  <- apply(R2_EFM0.99, 2, sd)
Mean_R2_EFM0.995 <- colMeans(R2_EFM0.995); sd_R2_EFM0.995  <- apply(R2_EFM0.995, 2, sd)
Mean_R2_EFM0.999 <- colMeans(R2_EFM0.999); sd_R2_EFM0.999  <- apply(R2_EFM0.999, 2, sd)

RMSE_QFM0.95  <- mean(recons_QFM0.95); SD_QFM0.95  <- sd(recons_QFM0.95)
RMSE_QFM0.99  <- mean(recons_QFM0.99); SD_QFM0.99  <- sd(recons_QFM0.99)
RMSE_QFM0.995 <- mean(recons_QFM0.995); SD_QFM0.995 <- sd(recons_QFM0.995)
RMSE_QFM0.999 <- mean(recons_QFM0.999); SD_QFM0.999 <- sd(recons_QFM0.999)

RMSE_EFM0.95  <- mean(recons_EFM0.95); SD_EFM0.95  <- sd(recons_EFM0.95)
RMSE_EFM0.99  <- mean(recons_EFM0.99); SD_EFM0.99  <- sd(recons_EFM0.99)
RMSE_EFM0.995 <- mean(recons_EFM0.995); SD_EFM0.995 <- sd(recons_EFM0.995)
RMSE_EFM0.999 <- mean(recons_EFM0.999); SD_EFM0.999 <- sd(recons_EFM0.999)

Mean_D_QFM0.95 <- mean(D_QFM0.95); SD_D_QFM0.95 <- sd(D_QFM0.95)
Mean_D_QFM0.99 <- mean(D_QFM0.99); SD_D_QFM0.99 <- sd(D_QFM0.99)
Mean_D_QFM0.995 <- mean(D_QFM0.995); SD_D_QFM0.995 <- sd(D_QFM0.995)
Mean_D_QFM0.999 <- mean(D_QFM0.999); SD_D_QFM0.999 <- sd(D_QFM0.999)

Mean_D_EFM0.95 <- mean(D_EFM0.95); SD_D_EFM0.95 <- sd(D_EFM0.95)
Mean_D_EFM0.99 <- mean(D_EFM0.99); SD_D_EFM0.99 <- sd(D_EFM0.99)
Mean_D_EFM0.995 <- mean(D_EFM0.995); SD_D_EFM0.995 <- sd(D_EFM0.995)
Mean_D_EFM0.999 <- mean(D_EFM0.999); SD_D_EFM0.999 <- sd(D_EFM0.999)

QFM0.95[[n]]  <- c(Mean_R2_QFM0.95,  RMSE_QFM0.95,  SD_QFM0.95, Mean_D_QFM0.95, SD_D_QFM0.95, sd_R2_QFM0.95)
QFM0.99[[n]]  <- c(Mean_R2_QFM0.99,  RMSE_QFM0.99,  SD_QFM0.99, Mean_D_QFM0.99, SD_D_QFM0.99,sd_R2_QFM0.99)
QFM0.995[[n]] <- c(Mean_R2_QFM0.995, RMSE_QFM0.995, SD_QFM0.995, Mean_D_QFM0.995, SD_D_QFM0.995,sd_R2_QFM0.995)
QFM0.999[[n]] <- c(Mean_R2_QFM0.999, RMSE_QFM0.999, SD_QFM0.999, Mean_D_QFM0.999, SD_D_QFM0.999, sd_R2_QFM0.999)

EFM0.95[[n]]  <- c(Mean_R2_EFM0.95,  RMSE_EFM0.95,  SD_EFM0.95, Mean_D_EFM0.95, SD_D_EFM0.95, sd_R2_EFM0.95)
EFM0.99[[n]]  <- c(Mean_R2_EFM0.99,  RMSE_EFM0.99,  SD_EFM0.99, Mean_D_EFM0.99, SD_D_EFM0.99, sd_R2_EFM0.99)
EFM0.995[[n]] <- c(Mean_R2_EFM0.995, RMSE_EFM0.995, SD_EFM0.995, Mean_D_EFM0.995, SD_D_EFM0.995, sd_R2_EFM0.995)
EFM0.999[[n]] <- c(Mean_R2_EFM0.999, RMSE_EFM0.999, SD_EFM0.999, Mean_D_EFM0.999, SD_D_EFM0.999, sd_R2_EFM0.999)

# Results
Table1_mixgauss = matrix(nrow=0, ncol=10) 
table = rbind(QFM0.95[[n]], EFM0.95[[n]], QFM0.99[[n]],EFM0.99[[n]], QFM0.995[[n]], EFM0.995[[n]], QFM0.999[[n]], EFM0.999[[n]])      


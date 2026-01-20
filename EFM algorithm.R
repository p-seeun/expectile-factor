library(expm)
rsquare <- function(a, B){
  lm<- summary(lm(a ~ B))
  return(lm$adj.r.squared)
}
efm_est = function(Y, r, tau, tol=1e-5){

Te = nrow(Y)
N = ncol(Y)

# 초기값 설정 
Lambda <- matrix(rnorm(N * r), nrow = N, ncol = r)
F <- matrix(rnorm(Te * r), nrow = Te, ncol = r)

# IRLS 시작
max_iter = 200
prev_loss = Inf
for (iter in 1:max_iter) {
  
  # Residual 계산 (Te x N)
  R <- Y- F %*% t(Lambda)

  # 가중치 계산
  W <- ifelse(R >= 0, tau, 1 - tau)
  
  # Lambda 업데이트 (F 고정)
  for (i in 1:N) {
    W_i <- diag(W[, i])  # w_{1i}, w_{2i}, \cdot, w_{Te i} 
    X <- F               # (Te x r)
    y <- Y[, i]          # (Te x 1)
    Lambda[i, ] <- solve(t(X) %*% W_i %*% X + 1e-6 * diag(r), t(X) %*% W_i %*% y)
  }
  
  # F 업데이트 (Lambda 고정)
  for (t_ in 1:Te) {
    W_t <- diag(W[t_, ])  # w_{t1}, w_{t2}, \cdot, w_{tN} 
    X <- Lambda           # (N x r)
    y <- Y[t_, ]          # (N x 1)
    F[t_, ] <- solve(t(X) %*% W_t %*% X + 1e-6 * diag(r), t(X) %*% W_t %*% y)
  }
  
  # Check convergence
  R_new <- Y - F %*% t(Lambda)
  loss <- mean(ifelse(R_new >= 0, tau, 1 - tau) * R_new^2)
  
  if (iter > 1 && abs(prev_loss - loss) < tol) {
    cat("Converged at iteration:", iter,"\n")
    break
  }
  prev_loss <- loss
}

# Normalization 
Fhat <- F
Lhat <- Lambda
sigmaF <- t(Fhat) %*% Fhat / Te
sigmaA <- t(Lhat) %*% Lhat / N

dum1 <- sqrtm(sigmaF) %*% sigmaA %*% sqrtm(sigmaF)
dum2 <- eigen(dum1)$vectors

R<- solve(sqrtm(sigmaF), dum2)
Fhat <- Fhat %*% R
Lhat <- Lhat %*% t(solve(R))

return(list(lmat = Lhat, fmat = Fhat))
}



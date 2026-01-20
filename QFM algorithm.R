library(Matrix)
library(expm)
library(quantreg)

loss.fun <-function(tau, vec){
  loss = sapply(vec, function(x){ return(ifelse( x>0, tau*x, (tau-1)*x)) }) 
  return( mean(loss) )
}

# QFM estimation by Chen et al.(2021)
qfm_est <- function(X, r, tol=1e-5, tau) {
  Te <- nrow(X)
  N <- ncol(X)
  
  F0 <- matrix(rnorm(Te * r), nrow = Te)
  
  obj <- 0
  newobj <- 1
  
  while (abs(newobj - obj) > tol) {
    lambda0 <- matrix(0, nrow=N, ncol=r)
    for (i in 1:N) {
      lambda0[i,]<-rq(X[,i] ~ F0-1, tau = tau)$coeff
    }
    obj <- loss.fun(tau, as.vector(X - F0 %*% t(lambda0)))
    
    F1 <- matrix(0,nrow = Te, ncol = r)
    for (j in 1:Te) {
      F1[j,] <- rq(X[j,] ~ lambda0-1, tau=tau)$coeff
    }
    F0<-F1  
    
    lambda1 <- matrix(nrow = N, ncol = r)
    for (i in 1:N) {
      lambda1[i,]<-rq(X[,i] ~ F1-1, tau = tau)$coeff
    }
    newobj <- (loss.fun(tau, as.vector(X - F1 %*% t(lambda1))))
  }
  
  Fhat <- F1
  Lhat <- lambda1
  
  sigmaF <- t(Fhat) %*% Fhat / Te
  sigmaA <- t(Lhat) %*% Lhat / N
  
  dum1 <- sqrtm(sigmaF) %*% sigmaA %*% sqrtm(sigmaF)
  dum2 <- eigen(dum1)$vectors
  
  R<- solve(sqrtm(sigmaF), dum2)
  Fhat <- Fhat %*% R
  Lhat <- Lhat %*% t(solve(R))
  
  return(list(fmat = Fhat, lmat = Lhat))
}


#### (아마 필요없을것)
loss.fun <-function(tau, vec){
  return(mean(sapply(vec,function(x){return(ifelse(x>0, tau*x, (1-tau)*abs(x)))})))
}

do_qfm <- function(X,r.vec,tau.vec){ #X should be centered #tau.vec needed 
  N = ncol(X)
  Te = nrow(X)
  temp <- svd(X)
  pred.quant.i = list()
  for(ind in 1:length(tau.vec)){
    #print(paste("###########", ind, "th QFM estimation###########", sep=""))
    tau <- tau.vec[ind] 
    r = ifelse(length(r.vec)==1,r.vec,r.vec[ind])
    L.hat = sqrt(N)*(temp$v)[,1:r] #initial
    F.hat = 1/N * X %*% L.hat      #initial
    
    # mean(abs(X-F.hat%*%t(L.hat)))
    L.hat.new <- matrix(0, nrow=N, ncol=r)
    F.hat.new <- matrix(0, nrow=Te, ncol=r)
    
    iter <- 1
    while(1){
      for(i in 1:N){
        L.hat.new[i,] <- rq(X[,i] ~ F.hat-1, tau=tau)$coeff
      }
      for(t in 1:nrow(F.hat)){
        F.hat.new[t,] <- rq(X[t,] ~ L.hat.new-1, tau=tau)$coeff
      }
      q1 = qr(F.hat.new)
      q2 = qr(qr.R(q1)%*%t(L.hat.new))
      F.hat.new = sqrt(Te) * qr.Q(q1) %*% qr.Q(q2)
      L.hat.new = t(qr.R(q2))/sqrt(Te)
      
      if(iter>1 ){
        if(abs(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector(X-F.hat %*% t(L.hat))))< 1e-3 ){
          #print(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector( X-F.hat %*% t(L.hat))))
          #print(iter)
          break
        }
      }
      iter <- iter+1
      F.hat <- F.hat.new
      L.hat <- L.hat.new
      if(iter==200){
        #print(paste("Max iteration at tau=", tau))
        break;
      }
    }
    pred.quant.i [[ind]] = list(F=F.hat.new, L=L.hat.new)
  }
  return(pred.quant.i)
}


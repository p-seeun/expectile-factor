# Forecasting function when factors are included
Frcst = function(Z.win, Y.win, h, pmax, Fact){
  Te = length(Y.win) 
  Y.st = matrix(0,nrow=Te-pmax-h+1, ncol=0)
  for(p in pmax:1){
    Y.st = cbind(Y.st, Y.win[p:(Te-h+p-pmax)])
  }
  BICs = c()
  # For p = 0 
  BICs[1] = BIC( lm(Z.win[(pmax+h):Te] ~ Fact[pmax:(Te-h),]) )
  # For p > 0
  for(p in 1:pmax){
    lm = lm(Z.win[(pmax+h):Te] ~ Y.st[, 1:p] + Fact[pmax:(Te-h),])
    BICs[p+1] = BIC(lm)
  }
  p.opt = which.min(BICs) - 1
  
  # If phat = 0
  if(p.opt==0){
    lm = lm(Z.win[(pmax+h):Te] ~ Fact[pmax:(Te-h),] )
    pred = lm$coefficients %*% c(1, Fact[Te,])
  }
  # If phat > 0
  else{
    lm = lm(Z.win[(pmax+h):Te] ~ Y.st[, 1:p.opt] + Fact[pmax:(Te-h),] )
    pred = lm$coefficients %*% c(1,Y.win[Te:(Te-p.opt+1)], Fact[Te,])
  }
  return(list(pred=pred, p.opt=p.opt))
}

#Forecasting function when only lags are included
Frcst_ar = function(Z.win, Y.win, h, pmax){
  Te = length(Y.win) 
  Y.st = matrix(0,nrow=Te-pmax-h+1, ncol=0)
  for(p in pmax:1){
    Y.st = cbind(Y.st, Y.win[p:(Te-h+p-pmax)])
  }
  BICs = c()
  # For p = 0
  BICs[1] = BIC( lm(Z.win[(pmax+h):Te] ~ rep(1,Te-pmax-h+1) ) )
  # For p > 0
  for(p in 1:pmax){
    lm = lm(Z.win[(pmax+h):Te] ~ Y.st[, 1:p] )
    BICs[p+1] = BIC(lm)
  }
  p.opt = which.min(BICs) - 1
  
  # If phat = 0
  if(p.opt==0){
    lm = lm(Z.win[(pmax+h):Te] ~ rep(1,Te-pmax-h+1) )
    pred = lm$coefficients[1]
  }
  # If phat > 0
  else{
    lm = lm(Z.win[(pmax+h):Te] ~ Y.st[, 1:p.opt] )
    pred = lm$coefficients %*% c(1,Y.win[Te:(Te-p.opt+1)])
  }
  return(list(pred=pred, p.opt=p.opt))
}


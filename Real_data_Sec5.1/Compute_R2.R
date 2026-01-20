source("../EFM algorithm.R")
source("../QFM algorithm.R")
library(quantmod)

#daily return
getSymbols("^GSPC", from = "2003-01-01", to = "2023-12-31")
spx <- na.omit(Cl(GSPC))
spxret <- na.omit(spx / lag(spx) - 1)*100

eval_dates <- seq(as.Date("2004-01-01"), as.Date("2023-12-31"), by = "month")

library(rugarch)

# Compute 5% VaR and 5% ES
fhs_garch <- function(ret_xts, eval_date,  alpha  = 0.1, n_sim  = 5000,  window = NULL) {
  set.seed(0)
  eval_date <- as.Date(eval_date)
  all_dates <- index(ret_xts)
  
  ref_date <- max(all_dates[all_dates <= eval_date])

  ## training data
  train_ret_xts <- ret_xts[all_dates <= ref_date]

  if (!is.null(window) && NROW(train_ret_xts) > window) {
    train_ret_xts <- tail(train_ret_xts, window)
  }
  r <- as.numeric(train_ret_xts)
  
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(1,1), include.mean = TRUE),
    distribution.model = "std"
  )
  
  fit <- ugarchfit(spec, data = r, solver = "hybrid")
  z_hist <- as.numeric(residuals(fit, standardize = TRUE))
  sigma_all <- as.numeric(sigma(fit))

  sigma_last <- tail(sigma_all, 1)
  mu_last   <- tail(as.numeric(fitted(fit)), 1)
  
  ## generate simulated return distribution
  sim_ret <- numeric(n_sim)
  
  for (k in 1:n_sim) {
    z_t   <- sample(z_hist, size = 1, replace = TRUE)
    r_new <- sigma_last *z_t + mu_last
    sim_ret[k] <- r_new
  }
  
  ## calculate VaR and ES 
  VaR <- as.numeric(quantile(sim_ret, probs = alpha, na.rm = TRUE))
  ES  <- mean(sim_ret[sim_ret <= VaR], na.rm = TRUE)
  
  out <- list(VaR = VaR, ES = ES)
  return(out)
}

var_es_monthly = data.frame(ES=rep(NA,240), VaR=rep(NA,240))
for(i in 1:length(eval_dates)){
  temp = fhs_garch(spxret, eval_date = eval_dates[i], alpha = 0.05, window = 150)
  var_es_monthly$VaR[i]= temp$VaR
  var_es_monthly$ES[i] = temp$ES
}
#############################
X = readRDS("CRSP(240months).rds")

# Estimate Mean Factor
svd = svd(scale(X, scale=FALSE))
N = nrow(X)
lambda = t(as.matrix(sqrt(N)*t(svd$v)[1:6,])) 
Factor = (scale(X, scale=FALSE) %*% lambda)/N 

#### Compute R^2

#QFM, EFM estimators are obtained from 'Get_EFM&QFM_estimators.R'
QFM = readRDS("QFM_estimators.rds")
EFM = readRDS("EFM_estimators.rds")

tau.vec= c(0.001, 0.005, seq(0.01,0.1,0.01))

dev.off()
####################
##### R^2 for VaR
####################

QFM_VaR = data.frame(Tot=rep(NA,length(tau.vec)), Covid=rep(NA,length(tau.vec)), Crisis=rep(NA,length(tau.vec)))
EFM_VaR = data.frame(Tot=rep(NA,length(tau.vec)), Covid=rep(NA,length(tau.vec)), Crisis=rep(NA,length(tau.vec)))

# 1) Total (2004.01-2023.12)
for(i in seq_along(tau.vec)){
  QFM_VaR$Tot[i]=round(rsquare(var_es_monthly$VaR, cbind(Factor,QFM[[i]]$fmat)),3)
}
for(i in seq_along(tau.vec)){
  EFM_VaR$Tot[i]=round(rsquare(var_es_monthly$VaR, cbind(Factor,EFM[[i]]$fmat)),3)
}
  
# 2) Covid (2020.01-2022.12)
period = which(eval_dates=="2020-01-01"):which(eval_dates=="2022-12-01")
for(i in seq_along(tau.vec)){
  QFM_VaR$Covid[i]=round(rsquare(var_es_monthly$VaR[period], cbind(Factor[period,],QFM[[i]]$fmat[period,])),3)
}
for(i in seq_along(tau.vec)){
  EFM_VaR$Covid[i]=round(rsquare(var_es_monthly$VaR[period], cbind(Factor[period,],EFM[[i]]$fmat[period,])),3)
}
  
# 3) Crisis (2007.07-2009.03)
period = which(eval_dates=="2007-07-01"):which(eval_dates=="2009-06-01")
for(i in seq_along(tau.vec)){
  QFM_VaR$Crisis[i]=round(rsquare(var_es_monthly$VaR[period], cbind(Factor[period,],QFM[[i]]$fmat[period,])),3)
}
for(i in seq_along(tau.vec)){
  EFM_VaR$Crisis[i]=round(rsquare(var_es_monthly$VaR[period], cbind(Factor[period,],EFM[[i]]$fmat[period,])),3)
}

####################
##### R^2 for ES
####################

QFM_ES = data.frame(Tot=rep(NA,length(tau.vec)), Covid=rep(NA,length(tau.vec)), Crisis=rep(NA,length(tau.vec)))
EFM_ES = data.frame(Tot=rep(NA,length(tau.vec)), Covid=rep(NA,length(tau.vec)), Crisis=rep(NA,length(tau.vec)))

# 1) Total (2004.01-2023.12)
for(i in seq_along(tau.vec)){
  QFM_ES$Tot[i]=round(rsquare(var_es_monthly$ES, cbind(Factor,QFM[[i]]$fmat)),3)
}
for(i in seq_along(tau.vec)){
  EFM_ES$Tot[i]=round(rsquare(var_es_monthly$ES, cbind(Factor,EFM[[i]]$fmat)),3)
}
  
# 2) Covid (2020.01-2022.12)
period = which(eval_dates=="2020-01-01"):which(eval_dates=="2022-12-01")
for(i in seq_along(tau.vec)){
  QFM_ES$Covid[i]=round(rsquare(var_es_monthly$ES[period], cbind(Factor[period,],QFM[[i]]$fmat[period,])),3)
}
for(i in seq_along(tau.vec)){
  EFM_ES$Covid[i]=round(rsquare(var_es_monthly$ES[period], cbind(Factor[period,],EFM[[i]]$fmat[period,])),3)
}
  
# 3) Crisis (2007.07-2009.03)
period = which(eval_dates=="2007-07-01"):which(eval_dates=="2009-06-01")
for(i in seq_along(tau.vec)){
  QFM_ES$Crisis[i]=round(rsquare(var_es_monthly$ES[period], cbind(Factor[period,],QFM[[i]]$fmat[period,])),3)
}
for(i in seq_along(tau.vec)){
  EFM_ES$Crisis[i]=round(rsquare(var_es_monthly$ES[period], cbind(Factor[period,],EFM[[i]]$fmat[period,])),3)
}

############################
#### Six plots on Figures 2 and 3
############################
par(mar = c(5, 4.7, 4, 2)) 
ylim_range <- range(EFM_VaR$Tot, QFM_VaR$Tot) +c(-0.1,0.1)
plot(tau.vec, EFM_VaR$Tot, type='l', ylim=ylim_range,col="red", ylab=expression(R^2), main="Total (VaR)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_VaR$Tot, col="blue")
points(tau.vec, EFM_VaR$Tot, col = "red", pch = 16)
points(tau.vec, QFM_VaR$Tot, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")

ylim_range <- range(EFM_ES$Tot, QFM_ES$Tot) +c(-0.1,0.1)
plot(tau.vec, EFM_ES$Tot, type='l', ylim=ylim_range,col="red", ylab=expression(R^2), main="Total (ES)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_ES$Tot, col="blue")
points(tau.vec, EFM_ES$Tot, col = "red", pch = 16)
points(tau.vec, QFM_ES$Tot, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")

ylim_range <- range(EFM_VaR$Covid, QFM_VaR$Covid) +c(-0.1,0.1)
plot(tau.vec, EFM_VaR$Covid, type='l',col="red", ylim = ylim_range, ylab=expression(R^2), main="COVID (VaR)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_VaR$Covid, col="blue")
points(tau.vec, EFM_VaR$Covid, col = "red", pch = 16)
points(tau.vec, QFM_VaR$Covid, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")

ylim_range <- range(EFM_ES$Covid, QFM_ES$Covid) +c(-0.1,0.1)
plot(tau.vec, EFM_ES$Covid, type='l',col="red", ylim = ylim_range, ylab=expression(R^2), main="COVID (ES)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_ES$Covid, col="blue")
points(tau.vec, EFM_ES$Covid, col = "red", pch = 16)
points(tau.vec, QFM_ES$Covid, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")

ylim_range <- range(EFM_VaR$Crisis, QFM_VaR$Crisis)+c(-0.1,0.1)
plot(tau.vec, EFM_VaR$Crisis, type='l',col="red", ylim = ylim_range,ylab=expression(R^2), main="Global crisis (VaR)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_VaR$Crisis, col="blue")
points(tau.vec, EFM_VaR$Crisis, col = "red", pch = 16)
points(tau.vec, QFM_VaR$Crisis, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")  

ylim_range <- range(EFM_ES$Crisis, QFM_ES$Crisis)+c(-0.1,0.1)
plot(tau.vec, EFM_ES$Crisis, type='l',col="red", ylim = ylim_range, ylab=expression(R^2), main="Global crisis (ES)", xlab="level", cex.lab=1.2)
lines(tau.vec, QFM_ES$Crisis, col="blue")
points(tau.vec, EFM_ES$Crisis, col = "red", pch = 16)
points(tau.vec, QFM_ES$Crisis, col = "blue", pch = 16)
legend("bottomright", legend=c("EFM","QFM"),
       col=c("red","blue"),
       pch=16, lty=1,
       bty="o")






source("../EFM algorithm.R")
source("../QFM algorithm.R")

#CRSP monthly data of 1382 stocks from Jan 2004 to Dec 2023.
X = readRDS("CRSP(240months).rds") 

## Estimate Quantile Factors
tau.vec= c(0.001, 0.005, seq(0.01,0.1,0.01))
r.vec = rep(1,12) #Factor numbers estimated from 'Facnum.R'

QFM=list()
for(i in 1:length(tau.vec)){
  set.seed(1)
  QFM[[i]] = qfm_est( scale(X, scale=FALSE), r=r.vec[i], tau=tau.vec[i] )
  print(i)
}
#saveRDS(QFM, "QFM_estimators.rds")


## Estimate Expectile Factors
tau.vec= c(0.001, 0.005, seq(0.01,0.1,0.01))
r.vec = c(1,4,4,3,3,3,3,2,2,2,4,4)

EFM=list()
for(i in 1:length(tau.vec)){
  set.seed(0)
  EFM[[i]] = efm_est( scale(X, scale=FALSE), r=r.vec[i], tau=tau.vec[i] )
  print(i)
}
#saveRDS(EFM, "EFM_estimators.rds")

source("../EFM algorithm.R")

X = readRDS("FRED-MD.rds") 

EFM1 = list()
for(st in 1:(nrow(X)-120)){
  set.seed(st)
  
  X.win = X[st:(st+119),]
  
  EFM1[[st]] = efm_est(X.win, r=5, tau=0.01, tol=1e-5)$fmat
  
  print(st)
}
#saveRDS(EFM1, "EFM_estimators/EFM1.rds")


EFM5 = list()
for(st in 1:(nrow(X)-120)){
  set.seed(st)

  X.win = X[st:(st+119),]

  EFM5[[st]] = efm_est(X.win, r=6, tau=0.05, tol=1e-5)$fmat
  
  print(st)
}
#saveRDS(EFM5, "EFM_estimators/EFM5.rds")

EFM50 = list()
for(st in 1:(nrow(X)-120)){
  set.seed(st)

  X.win = X[st:(st+119),]
  
  EFM50[[st]] = efm_est(X.win, r=5, tau=0.5, tol=1e-5)$fmat
  
  print(st)
}
#saveRDS(EFM50, "EFM_estimators/EFM50.rds")

EFM95 = list()
for(st in 1:(nrow(X)-120)){
  set.seed(st)

  X.win = X[st:(st+119),]
  
  EFM95[[st]] = efm_est(X.win, r=7, tau=0.95, tol=1e-5)$fmat
  
  print(st)
}
#saveRDS(EFM95, "EFM_estimators/EFM95.rds")


EFM99 = list()
for(st in 1:(nrow(X)-120)){
  set.seed(st)
  
  X.win = X[st:(st+119),]
  
  EFM99[[st]] = efm_est(X.win, r=6, tau=0.99, tol=1e-5)$fmat
  
  print(st)
}
#saveRDS(EFM99, "EFM_estimators/EFM99.rds")




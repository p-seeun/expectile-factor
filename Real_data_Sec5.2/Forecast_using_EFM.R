source("../EFM algorithm.R")
source("Forecast functions.R")

X = readRDS("FRED-MD.rds") 

library(readr)
dum <- read_csv("2024-06.csv")
dum$INDPRO
which(dum$sasdate == "1/1/1990") #374
which(dum$sasdate == "12/1/2019") #733

Y = log(dum$INDPRO[374:733]/dum$INDPRO[(374-1):(733-1)])

EFM1 = readRDS("EFM_estimators/EFM1.rds")
EFM5 = readRDS("EFM_estimators/EFM5.rds")
EFM50 = readRDS("EFM_estimators/EFM50.rds")
EFM95 = readRDS("EFM_estimators/EFM95.rds")
EFM99 = readRDS("EFM_estimators/EFM99.rds")

pmax = 5

for(h in 1:6){
  lags = log(dum$INDPRO[374:733] / dum$INDPRO[(374-h):(733-h)])
  Y.true = exp(lags[(120+h):240]) * dum$INDPRO[(733-120):(733-h)]
  
  # AR
  pred = c()
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]
    pred[st] = Frcst_ar(Z.win, Y.win, h, pmax)$pred
  }
  pred = exp(pred) * dum$INDPRO[(733-120):(733-h)]
  ar = sqrt(mean((pred-Y.true)^2)) 
  
  pred1=c(); pred5=c(); pred50=c(); pred95=c(); pred99=c(); 
  
  for(st in 1:(120-h+1)){
    Z.win = lags[st:(st+119)]
    Y.win = Y[st:(st+119)]

    res1 = EFM1[[st]]
    res5 = EFM5[[st]]
    res50 = EFM50[[st]]
    res95 = EFM95[[st]]
    res99 = EFM99[[st]]
    
    pred1[st] = Frcst(Z.win, Y.win, h, pmax, Fact=res1)$pred
    pred5[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(res5))$pred
    pred50[st] = Frcst(Z.win, Y.win, h, pmax, Fact=res50)$pred
    pred95[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(res95))$pred
    pred99[st] = Frcst(Z.win, Y.win, h, pmax, Fact=cbind(res99))$pred
  }
 
  pred1 = exp(pred1) * dum$INDPRO[(733-120):(733-h)]
  pred5 = exp(pred5) * dum$INDPRO[(733-120):(733-h)]
  pred50 = exp(pred50) * dum$INDPRO[(733-120):(733-h)]
  pred95 = exp(pred95) * dum$INDPRO[(733-120):(733-h)]
  pred99 = exp(pred99) * dum$INDPRO[(733-120):(733-h)]
  
  efm1 = sqrt(mean((pred1-Y.true)^2))
  efm5 = sqrt(mean((pred5-Y.true)^2))
  efm50 = sqrt(mean((pred50-Y.true)^2))
  efm95 = sqrt(mean((pred95-Y.true)^2))
  efm99 = sqrt(mean((pred99-Y.true)^2))
  
  print(h)
  lst = list(ar1=ar, ar=(ar/ar)^2, efm1=(efm1/ar)^2, efm5=(efm5/ar)^2, efm50=(efm50/ar)^2, efm95=(efm95/ar)^2, efm99=(efm99/ar)^2)  
  print(lapply(lst, function(x) round(x, 5)))
}



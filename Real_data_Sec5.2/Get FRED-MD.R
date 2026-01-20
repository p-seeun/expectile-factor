#devtools::install_github("cykbennie/fbi")
library(fbi)
DATA = fredmd(file="2024-06.csv", transform=TRUE) #transformed to be stationary
dates = DATA[,1]
data = DATA[,-1]

dates <- dates[373:732]  # 1990-01 to 2019-12
data <- data[373:732, ]

# Remove columns with all missing values
X_raw <- data[, colSums(is.na(data)) == 0]  
dim(X_raw) # 360*125

X = scale(X_raw)
#saveRDS(X, "FRED-MD.rds")




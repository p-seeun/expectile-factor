# expectile-factor
This repository provides the source code for simulations and real data analyses included in the paper 'Tail-oriented Factor Models Using Expectiles' is available. For the real data analyses, two economic datasets are used: the CRSP monthly stock return data and the FRED-MD data.

## Overview
- `EFM algorithm.R`: It is a code for functions that estimate expectile factor model. 
- `QFM algorithm.R`: It is a code that implements algorithms of quantile factor model proposed by Chen et al.(2021).
- Simulations
  - This folder includes two R code files that generate the simulation results for Section 4.1 and Section 4.2, respectively.
- Real_data_Sec5.1
  - This folder includes codes that are needed to analyze CRSP monthly stock return data in Section 5.1. Unfortunately, the dataset is not freely available, but many institutions have access to WRDS, which provides CRSP data. The code assumes that the monthly stock return data of 1,382 stocks, from January 2004 to December 2023, are saved as 'CRSP(240months).rds'.
  - `Factor_num.R` : It estimates the number of factors of EFM and QFM across different levels.
  - `Get_EFM&QFM_estimators.R` : This code obtains the estimators of EFM and QFM applied to the dataset across different levels. The resulting estimators are saved as `EFM_estimators.rds` and `QFM_estimators.rds`.
  - `Compute_R2.R` : It computes the adjusted $R^2$ values and generates all results presented in Figures 2 and 3.
- Real_data_Sec5.2
  - This folder includes codes for empirical applications using FRED-MD dataset in Section 5.2.
  - `Get FRED-MD.R` : It processes the panel data of FRED-MD data from the raw data `2024-06.csv`, and the results are saved as `FRED-MD.rds`.
  - `Factor_num.R` : It estimates the number of factors of EFM across different levels.
  - `Get_EFM_estimators.R` : The factors from EFM of different levels on each rolling-window data are estimated and saved as .rds files in the folder `EFM_estimators`.
  - `Forecast_using_EFM.R` : The main code that forecasts the industrial production index and derives results in Table 6.
  - `Forecast functions.R` : It includes functions needed in `Forecast_using_EFM.R`.
  - `Draw_heatmap.R` : This code draws heatmaps presented in Figure 4.
  

#clear environment
rm(list=ls(all=T)) 
library(tictoc)

#path
setwd("~/Dropbox/Research/TIBs/Rcodes_RFS") # "/Users/pchikhan/Dropbox/TIBs/Rcodes_v8_CLRP" 

#Optimization of parameters
indic_estim       <-0 # 1 if you want to launch an estimation
indic_plots_paper <-1 # 1 if you want to run some plots
indic_tables_paper<-0 # 1 if you want to run tables

#For scripts using parallel computing:
nb.cores <- 8 #number of cores you want to dedicate to parallel computations

if(indic_estim==1){
  print("Optimize the model: it might take few minutes...")
  source("./estimations/load_ini_model.R")
  source("./estimations/plots_check.R")
  source("./estimations/optimiz.R")
  source("./estimations/plots_check.R")
}else{
  print("Load optimized model: it might take few seconds...")
  source("./estimations/load_ini_model.R")
  source("./estimations/optimiz.R")
  source("./estimations/plots_check.R")
}

#Updating Plots and Tables

##Plots
#* Description
#* 1 tibs and swaps
#* 2 digital option
#* 3 SCC-RP
#* 4 pdf temperatures
#* 5 climate beta
#* 6 disasters simulations
#* 7 pdf sea level
#* 8 mu (comparison with DICE)
#* 9 Radiative forcing approximation
#*10 constant maturity - ZCB
#*11 cut Climate Premium

plots <-c(2) #Choose here the plots you would like to run, from 1:11


##Tables
#* Description
#* -Estimated parameters
#* -Target vs model-implied moments
#* -List of parameters
#* -Initial values of state vector


if(indic_plots_paper==1){
  source("outputs/plots_paper.R")
}

if(indic_tables_paper==1){
  source("outputs/tables_paper.R")
}
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


library(parallel)
library(doParallel)
library(mgcv)
library("colorspace")
library("broman")

# Maturity:
H <-model_sol$horiz.2100                                                        #maturity 2100
EV<-EV.fct(model_sol,H)
#Pricing
omega_ZCB  <- model_sol$omega_ZCB
omega_T.at <- model_sol$omega_T.at

x <- exp(seq(-20,40,length.out = 2000))                                         #grid for Proposition 8 (Fourier)
a <- omega_T.at                                                                 #Options related to T_at in X
b <- 2                                                                          #Options pay when a'X < b

# Miscellaneous
n.date<-model_sol$Tmax

# For scripts using parallel computing:
number.of.cores <- nb.cores

# Vector of temperatures for strike K, Digital Option:
K <- c(2.5,3,4)

# Quantiles considered for pdf charts:
vector.of.CI <- c(0.5,0.8,0.9,0.95)

# Colors:
P.col.line <- brocolors("crayons")["Teal Blue"]
P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
Q.col.line <- brocolors("crayons")["Mango Tango"]
Q.col<-adjustcolor( Q.col.line, alpha.f = 0.15)


print("Indicated time for charts: 8 cores.")
#********************************1*********************************************#
#TIBS/SWAPS
if(is.element(1,plots)){
  print("Producing plot 1: Swaps and TIBs - 1 second...")
  source("outputs/make_figures/make_figure_TS.R")
}

#********************************2*********************************************#
#Digital Option
if(is.element(2,plots)){
  print("Producing plot 2: Digital Options - 1 minute...")
  source("outputs/make_figures/make_figure_options.R")
}

#********************************3*********************************************#
#SCC-RP
if(is.element(3,plots)){
  print("Producing plot 3: Contour plot - 30 seconds...")
  print("Parallel computations - 8 cores")
  source("outputs/make_figures/make_figure_contourplots.R")
}

#********************************4*********************************************#
#Temperatures pdf + RCP
if(is.element(4,plots)){
  print("Producing plot 4: Atm. temperature PDF - 2 minutes...")
  print("Parallel computations - 8 cores")
  source("outputs/make_figures/make_figure_Tpdf.R")
}

#********************************5*********************************************#
#Climate beta
if(is.element(5,plots)){
  print("Producing plot 5: Climate beta - 1 minute...")
  print("Parallel computations - 8 cores")
  source("outputs/make_figures/make_figure_ClimateBeta.R")
}

#********************************6*********************************************#
#Disasters Simulations
if(is.element(6,plots)){
  print("Producing plot 6: Damages simulations - 1 second...")
  source("outputs/make_figures/make_figure_Dsimul.R")
}

#********************************7*********************************************#
#Global Sea Level pdf + RCP
if(is.element(7,plots)){
  print("Producing plot 7: Sea level rise PDF - 4 minutes...")
  print("Parallel computations - 8 cores")
  source("outputs/make_figures/make_figure_Hpdf.R")
}

#********************************8*********************************************#
#Mitigation vs DICE2016
if(is.element(8,plots)){
  print("Producing plot 8: Mitigation rate comparison - 1 second...")
  source("outputs/make_figures/make_figure_mu.R")
}

#********************************9*********************************************#
#Radiative Forcings Approximation
if(is.element(9,plots)){
  print("Producing plot 9: RF approximation - 1 second...")
  source("outputs/make_figures/make_figure_RFapprox.R")
}

#*******************************10*********************************************#
#Constant maturity for ZCB
if(is.element(10,plots)){
  print("Producing plot 10: Constant maturity rate ZCB - 3 seconds...")
  source("outputs/make_figures/make_figure_ConstantMaturityZCB.R")
}

#*******************************11*********************************************#
#Cut in Climate Premium
if(is.element(11,plots)){
  print("Producing plot 11: Climate Premium, cut in contour plot - 15 seconds...")
  print("Parallel computations - 8 cores")
  source("outputs/make_figures/make_figure_cut_CP_muD.R")
}

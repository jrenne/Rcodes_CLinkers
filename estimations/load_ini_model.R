#Initial model - solves model with initial parameters
source("./procedure/Functions_General_v2.R")
#library
library(optimx)
library(MASS)
library(expm)
##################
#Data information#
##################
#*X = delc, ytilde, E, Eind, F, Mat, Mup, Mlo, Tat, Tlo, CumD, CumE, Cumdc, H,
#*W

##################
#Data preparation#
##################
#state vector                                                                   #CHANGE
n.Z    <-14 
n.eta  <-3
n.W    <-n.eta+2 
Tmax   <-100
#eta
Phi.prep     <-matrix(0,n.eta,n.eta)
Phi.prep[2,2]<-0.95
Phi.prep[1,1]<-0
indic_stocks <-0
passive.mu   <-0
####################
#List of parameters#
####################
param<-list(
  A_bar      = 0.4,                                                             
  sigma_a    = 0.03,
  Phi        = Phi.prep,
  gamma      = 7,
  delta      = (1 - .015)^5,                                                    
  m0         = 1168/607,                                                        #RCP4.5+6.0+CDICE                                                      
  ell0.N     = 0.06,
  ell1.N     = matrix(0,n.Z+n.W,1),
  rho.N      = 0.01,
  ell0.D     = 0,                                                  
  ell1.D     = matrix(0,n.Z+n.W,1),                                                  
  mu_d       = 0.07,
  mu_n       = 6,
  eps_0      = 2.6,                                                             #DICE2016,eland0
  rho        = 0.115,                                                           #DICE2016,deland
  mateq      = 607,                                                             #CDICE
  mueq       = 600,                                                             #CDICE
  mleq       = 1772,                                                            #CDICE
  phi_0      = 0.5,                                                             #DICE2016,fex0
  phi_1      = 1,                                                               #DICE2016,fex1
  m_pi       = 607,                                                             #CDICE
  xi_1       = 0.685,                                                           #CDICE,c1
  xi_2       = 0.500,                                                           #CDICE,c3
  xi_3       = 0.034,                                                           #CDICE,c4
  lambda     = 3.45/3.25,                                                       #DICE2016,noname=fco22x/t2xco2
  tau        = 3.45,                                                            #DICE2016,fco22x
  delsigma   =-0.001,                                                           #DICE2016,dsig
  e0         = 35.85,                                                           #DICE2016
  q0         = 105.5,                                                           #DICE2016
  mu0        = 0.03,                                                            #DICE2016
  sigma0     = 35.85/(105.5*(1-0.03)),                                          #DICE2016
  gsigma1    =-0.0152,                                                          #DICE2016
  varphi_11  = 1-0.053,                                                         #CDICE,b11
  varphi_12  = 0.053,                                                           #CDICE,b12
  varphi_21  = 0.053*607/600,                                                   #CDICE,b21
  varphi_22  = 1-0.053*607/600-0.0042,                                          #CDICE,b22
  varphi_23  = 0.0042,                                                          #CDICE,b23
  varphi_32  = 0.0042*600/1772,                                                 #CDICE,b32
  varphi_33  = 1-0.0042*600/1772,                                               #CDICE,b33
  gback      = 0.025,                                                           #DICE2016
  pback      = 4000,
  theta2     = 2.6,                                                             #DICE2016,expcost2
  sigma_eta_f= 0.25, 
  weights_Tg = c(0.7,0.3),                                                      #Unused
  a_sat      = 0.003,                                                           #Rahmstorf&Vermeer(2009)
  T_0s       = 0.325,                                                           #Rahmstorf&Vermeer(2009)
  b_sat      = 0.025,                                                           #Rahmstorf&Vermeer(2009)
  delta_K    = 0.27,                                                            #Hall&Jones(1999)
  ini_delc   = 0,
  ini_tildey = 0,
  ini_E      = 38.45,
  ini_Eind   = 35.85,
  ini_F      = 2,
  ini_Mat    = 847,                                                             #Chris Smith
  ini_Mup    = 822,                                                             #CDICE
  ini_Mlo    = 1810,                                                            #CDICE
  ini_Tat    = 1.11,                                                            #Chris Smith
  ini_Tlo    = 0.25,                                                            #Chris Smith
  ini_CumD   = 0,
  ini_CumE   = 38.45,
  ini_Cumdelc= 0,
  ini_H      = 0.13,
  tol.GN     = 10^(-6),
  eps.GN     = 10^(-5),
  sigma_H    = 0,
  ell0.div   = 0,
  ell1.div   = matrix(0,n.Z+n.W,1),
  mu_div     = 0,
  chi.div    = 2,
  b_sk       = 0.3/80,                                                          #Hinkel et al.(2014)
  c_sat      = -0.01
  )
#CDICE: Folini et al.(2021)
#DICE2016: Nordhaus(2017)

#Determine dependence of shocks D and N
param$ell1.N[9]   <- 0.08
param$ell1.D[9]   <- 0.15
param$ell1.div[14]<- 0.10

#matrix of ell1.N.tilde for gamma0
ell1.N.tilde             <-param$ell1.N
ell1.N.tilde[n.Z+n.eta+2]<-param$rho.N/param$mu_n

#Optimization of parameters
#upper and lower bounds
#list with min/max
bounds<-list(
  min_A_bar    =  exp(.005*5)/param$delta-(1-param$delta_K),
  max_A_bar    =  exp(.05*5)/param$delta-(1-param$delta_K),
  min_sigma_a  =  .0000000001*(param$A_bar+1-param$delta_K),
  max_sigma_a  =  .10*(param$A_bar+1-param$delta_K),
  min_mu_d     =  0,
  max_mu_d     = .3,
  max_mu_n     = 30,
  max_ell0.N   = .15,
  min_ell1.N   = .05,
  max_ell1.N   = .15,
  max_rho.N    = .9,
  max_ell1.D   = .25,
  min_pback    = 400,
  max_pback    = 4500,
  min_a_sat    = 0.0005,
  max_a_sat    = 0.004,
  min_sigma_H  = 0,
  max_sigma_H  = 0.2
)

# These are the moments the model should replicate (as well as possible):
target_vector <- c(
  2.9,                                                                          #1. E(temperature) in 2100
  0.32,                                                                         #2. Standard deviation temp 2100
  0.20,                                                                         #3. Permafrost temp change 2100
  60,                                                                           #4. Permafrost Cum_E 2100
  -0.125,                                                                       #5. Slope (Cum_D, Temp)
  1,                                                                            #6. LT rate target
  0.45,                                                                         #7. E(SL) in 2100
  0.1,                                                                          #8. Standard deviation sl 2100
  4.9,                                                                          #9. E(F) in 2100
  5.2,                                                                          #10.E(F) in 2200
  5.2,                                                                          #11.E(F) in 2300
  3.5,                                                                          #12.E(temperature) in 2200
  3.5,                                                                          #13.E(temperature) in 2300
  0,                                                                            #14.LT rate above 0
  1.8,                                                                          #15.E(temperature) in 2040
  2.1                                                                           #16.E(temperature) in 2050
)
# Define weights (relative importance of the different targets' components):
weights    <- 1/unlist(target_vector)^2
weights[3] <- weights[4]

weights[5] <- 40*weights[5]
weights[6] <- 40*weights[6]

weights[1]    <- 50
weights[15:16]<- 100

weights[9:14] <- 50

#number of periods in each t
tstep<-5

#dates
vec_date<-seq(2020,by=tstep,len=Tmax-1)

#Initial consumption
c0<-299                                                                         #In trillions, consumption over 5y

#Optimization of mu                                                             #Sensible to initial conditions
# theta0  <-list(c(log(0.17),-1/21*log(0.17)),                                    #DICE
#                c(20,-1/21*log(0.17)),                                           #No action   (mu=0 for all t)
#                c(5,-1/21*log(0.17)))                                            #100% action (mu=1 for all t)

theta0  <-list(c(log(0.17),-1/21*log(0.17)))
MAXIT   <-200 # 500 #1000

#Pricing
omega_ZCB     <-matrix(0,n.Z+n.W,1)
omega_T.at    <-omega_ZCB
omega_T.at[9] <-1
omega_Div     <-omega_ZCB
#omega_Div[15] <-1

#Dividends & Stocks
mu_div.0 <-matrix(0,nrow=Tmax)
kap0.om0 <-matrix(0,nrow=Tmax) #kappa0
mu_r0.om0<-matrix(0,nrow=Tmax) #mu_r0
kap1.a0  <-matrix(0,nrow=Tmax) #kappa1
mu_pd.a0 <-rep(list(matrix(0,nrow=n.Z+n.W+1)),Tmax) #mu_pd
mu_r1.a1 <-rep(list(matrix(0,nrow=n.Z+n.W)),Tmax) #mu_r1
ini_pd   <-0

#Separating the TP
mu_c     <- matrix(0,nrow=Tmax)
indic_tp <- 0


#Solving the SDF
ini_matx<-list()
inf_matx<-list()
model<-list("parameters"=param,"vec_date"=vec_date,"c0"=c0,"tstep"=tstep,
            "MAXIT"=MAXIT,"n.eta"=n.eta,"n.W"=n.W,"n.Z"=n.Z,"Tmax"=Tmax,
            "ell1.N.tilde"=ell1.N.tilde,"theta0"=theta0,
            "bounds"=bounds,"target_vector"=target_vector,"weights"=weights,
            "omega_ZCB"=omega_ZCB,"omega_T.at"=omega_T.at,"omega_Div"=omega_Div,
            "mu_div.0"=mu_div.0,"kap0.om0"=kap0.om0,"mu_r0.om0"=mu_r0.om0,
            "kap1.a0"=kap1.a0,"mu_pd.a0"=mu_pd.a0,"mu_r1.a1"=mu_r1.a1,
            "ini_pd"=ini_pd,"ini_matx"=ini_matx,"inf_matx"=inf_matx,
            "mu_c"=mu_c)

remove(ell1.N.tilde,Phi.prep)
tic()
model_sol<-model_solve(model,theta0)
#Define the horizon for optimization
horiz<-(2100-model_sol$vec_date[1])/model_sol$tstep

model[["horiz.2100"]]<-horiz
toc()
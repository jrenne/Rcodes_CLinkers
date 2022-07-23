#Functions for Calibration
#Optimization and loss function

##Model from G_code
# Create function that converts a model to a vector of parameters:
Model2Param <- function(model){
  n.eta <- model$n.eta
  param <- model$parameters
  n.Z   <- model$n.Z
  n.W   <- model$n.W
  bounds<- model$bounds
  
  param_vector <- rep(0,68+n.eta+3*(n.Z+n.W))                                   #CHANGE nb param
  x<-(param$A_bar-bounds$min_A_bar)/(bounds$max_A_bar-bounds$min_A_bar)
  param_vector[1] <-log(x/(1-x))
  
  x<-(param$sigma_a-bounds$min_sigma_a)/(bounds$max_sigma_a-bounds$min_sigma_a)
  param_vector[2] <-log(x/(1-x))
  
  param_vector[3:(2+n.eta)]  <-log(diag(param$Phi)/(1 - diag(param$Phi)))
  param_vector[2+n.eta+1]    <-log(param$gamma)
  param_vector[2+n.eta+2]    <-log(param$delta)
  param_vector[2+n.eta+3]    <-log(param$m0)
  
  x <-param$ell0.N/bounds$max_ell0.N
  param_vector[2+n.eta+4]    <-log(x/(1-x))
  
  x<-(param$ell1.N-matrix(c(rep(0,8),bounds$min_ell1.N,rep(0,n.Z+n.W-9)),ncol=1))/
    (bounds$max_ell1.N-bounds$min_ell1.N)
  param_vector[(2+n.eta+5):(2+n.eta+4+n.Z+n.W)]<-matrix(log(x/(1-x)),
                                                    ncol=1)
  
  x <-param$rho.N/bounds$max_rho.N
  param_vector[6+n.eta+n.Z+n.W+1]<-log(x/(1-x))
  param_vector[6+n.eta+n.Z+n.W+2]<-log(param$ell0.D)
  
  x<-param$ell1.D/bounds$max_ell1.D
  param_vector[(6+n.eta+n.Z+n.W+3):(6+n.eta+2*(n.Z+n.W)+2)]<-matrix(log(x/(1-x)),
                                                            ncol=1)
  
  x<-(param$mu_d-bounds$min_mu_d)/(bounds$max_mu_d-bounds$min_mu_d)
  param_vector[8+n.eta+2*(n.Z+n.W)+1] <-log(x/(1-x))
  
  x<- param$mu_n/bounds$max_mu_n
  param_vector[8+n.eta+2*(n.Z+n.W)+2] <-log(x/(1-x))
  
  param_vector[8+n.eta+2*(n.Z+n.W)+3] <-log(param$eps_0)
  param_vector[8+n.eta+2*(n.Z+n.W)+4] <-log(param$rho)
  param_vector[8+n.eta+2*(n.Z+n.W)+5] <-log(param$mateq)
  param_vector[8+n.eta+2*(n.Z+n.W)+6] <-log(param$mueq)
  param_vector[8+n.eta+2*(n.Z+n.W)+7] <-log(param$mleq)
  param_vector[8+n.eta+2*(n.Z+n.W)+8] <-log(param$phi_0)
  param_vector[8+n.eta+2*(n.Z+n.W)+9] <-log(param$phi_1)
  param_vector[8+n.eta+2*(n.Z+n.W)+10]<-log(param$m_pi)
  param_vector[8+n.eta+2*(n.Z+n.W)+11]<-log(param$xi_1)
  param_vector[8+n.eta+2*(n.Z+n.W)+12]<-log(param$xi_2)
  param_vector[8+n.eta+2*(n.Z+n.W)+13]<-log(param$xi_3)
  param_vector[8+n.eta+2*(n.Z+n.W)+14]<-log(param$lambda)
  param_vector[8+n.eta+2*(n.Z+n.W)+15]<-log(param$tau)
  param_vector[8+n.eta+2*(n.Z+n.W)+16]<-log(-param$delsigma)
  param_vector[8+n.eta+2*(n.Z+n.W)+17]<-log(param$e0)
  param_vector[8+n.eta+2*(n.Z+n.W)+18]<-log(param$q0)
  param_vector[8+n.eta+2*(n.Z+n.W)+19]<-log(param$mu0)
  param_vector[8+n.eta+2*(n.Z+n.W)+20]<-log(param$sigma0)
  param_vector[8+n.eta+2*(n.Z+n.W)+21]<-log(-param$gsigma1)
  param_vector[8+n.eta+2*(n.Z+n.W)+22]<-log(param$varphi_11)
  param_vector[8+n.eta+2*(n.Z+n.W)+23]<-log(param$varphi_12)
  param_vector[8+n.eta+2*(n.Z+n.W)+24]<-log(param$varphi_21)
  param_vector[8+n.eta+2*(n.Z+n.W)+25]<-log(param$varphi_22)
  param_vector[8+n.eta+2*(n.Z+n.W)+26]<-log(param$varphi_23)
  param_vector[8+n.eta+2*(n.Z+n.W)+27]<-log(param$varphi_32)
  param_vector[8+n.eta+2*(n.Z+n.W)+28]<-log(param$varphi_33)
  param_vector[8+n.eta+2*(n.Z+n.W)+29]<-log(param$gback)
  
  x<-(param$pback-bounds$min_pback)/(bounds$max_pback-bounds$min_pback)
  param_vector[8+n.eta+2*(n.Z+n.W)+30]<-log(x/(1-x))
  param_vector[8+n.eta+2*(n.Z+n.W)+31]<-log(param$theta2)
  param_vector[8+n.eta+2*(n.Z+n.W)+32]<-log(param$sigma_eta_f)
  
  param_vector[(41+n.eta+2*(n.Z+n.W)):
                 (42+n.eta+2*(n.Z+n.W))]<-log(param$weights_Tg)
  
  x<-(param$a_sat-bounds$min_a_sat)/(bounds$max_a_sat-bounds$min_a_sat)
  param_vector[43+n.eta+2*(n.Z+n.W)]  <-log(x/(1-x))
  
  param_vector[44+n.eta+2*(n.Z+n.W)]  <-log(param$T_0s)
  param_vector[45+n.eta+2*(n.Z+n.W)]  <-log(param$b_sat)
  param_vector[46+n.eta+2*(n.Z+n.W)]  <-log(param$delta_K)
  param_vector[47+n.eta+2*(n.Z+n.W)]  <-log(param$ini_delc)
  param_vector[48+n.eta+2*(n.Z+n.W)]  <-log(param$ini_tildey)
  param_vector[49+n.eta+2*(n.Z+n.W)]  <-log(param$ini_E)
  param_vector[50+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Eind)
  param_vector[51+n.eta+2*(n.Z+n.W)]  <-log(param$ini_F)
  param_vector[52+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Mat)
  param_vector[53+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Mup)
  param_vector[54+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Mlo)
  param_vector[55+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Tat)
  param_vector[56+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Tlo)
  param_vector[57+n.eta+2*(n.Z+n.W)]  <-log(param$ini_CumD)
  param_vector[58+n.eta+2*(n.Z+n.W)]  <-log(param$ini_CumE)
  param_vector[59+n.eta+2*(n.Z+n.W)]  <-log(param$ini_Cumdelc)
  param_vector[60+n.eta+2*(n.Z+n.W)]  <-log(param$ini_H)
  param_vector[61+n.eta+2*(n.Z+n.W)]  <-log(param$tol.GN)
  param_vector[62+n.eta+2*(n.Z+n.W)]  <-log(param$eps.GN)
  
  x<-(param$sigma_H-bounds$min_sigma_H)/(bounds$max_sigma_H-bounds$min_sigma_H)
  param_vector[63+n.eta+2*(n.Z+n.W)]    <-log(x/(1-x))
  param_vector[64+n.eta+2*(n.Z+n.W)]    <-log(param$ell0.div)
  param_vector[(64+n.eta+2*(n.Z+n.W)+1):
                 (64+n.eta+3*(n.Z+n.W))]<-log(param$ell1.div)
  param_vector[65+n.eta+3*(n.Z+n.W)]    <-log(param$mu_div)
  param_vector[66+n.eta+3*(n.Z+n.W)]    <-log(param$chi.div)
  param_vector[67+n.eta+3*(n.Z+n.W)]    <-log(param$b_sk)
  param_vector[68+n.eta+3*(n.Z+n.W)]    <-log(-param$c_sat)
  return(param_vector)
}
# Inverse function: converts param to model:
Param2Model <- function(param_vector,model){
  new_model <- model
  n.eta     <- model$n.eta 
  n.Z       <- model$n.Z
  n.W       <- model$n.W
  bounds    <- model$bounds
  
  
  ell0.N     = bounds$max_ell0.N*exp(param_vector[2+n.eta+4])/
    (1+exp(param_vector[2+n.eta+4]))
  ell1.N     = matrix(c(rep(0,8),bounds$min_ell1.N,rep(0,n.Z+n.W-9)),ncol=1)+
    matrix((bounds$max_ell1.N-bounds$min_ell1.N)*
             exp(param_vector[(2+n.eta+5):(2+n.eta+4+n.Z+n.W)])/
             (1+exp(param_vector[(2+n.eta+5):(2+n.eta+4+n.Z+n.W)])),
           ncol=1)

  rho.N      = bounds$max_rho.N*exp(param_vector[6+n.eta+n.Z+n.W+1])/
    (1+exp(param_vector[6+n.eta+n.Z+n.W+1]))
  ell0.D     = exp(param_vector[6+n.eta+n.Z+n.W+2])
  ell1.D     = matrix(bounds$max_ell1.D*
                        exp(param_vector[(6+n.eta+n.Z+n.W+3):
                                           (6+n.eta+2*(n.Z+n.W)+2)])/
                        (1+exp(param_vector[(6+n.eta+n.Z+n.W+3):
                                              (6+n.eta+2*(n.Z+n.W)+2)])),
                      ncol=1)
  mu_d       = bounds$min_mu_d+
    (bounds$max_mu_d-bounds$min_mu_d)*
    exp(param_vector[8+n.eta+2*(n.Z+n.W)+1])/
    (1+exp(param_vector[8+n.eta+2*(n.Z+n.W)+1]))
  mu_n       = bounds$max_mu_n*exp(param_vector[8+n.eta+2*(n.Z+n.W)+2])/
    (1+exp(param_vector[8+n.eta+2*(n.Z+n.W)+2]))
  
  ell1.N.tilde             <-ell1.N
  ell1.N.tilde[n.Z+n.eta+2]<-rho.N/mu_n
  
  A_bar      = bounds$min_A_bar+
    (bounds$max_A_bar-bounds$min_A_bar)*
    exp(param_vector[1])/(1+exp(param_vector[1]))
  
  sigma_a    = bounds$min_sigma_a+
    (bounds$max_sigma_a-bounds$min_sigma_a)*
    exp(param_vector[2])/(1+exp(param_vector[2]))
  
  pback      = bounds$min_pback+
    (bounds$max_pback-bounds$min_pback)*
    exp(param_vector[38+n.eta+2*(n.Z+n.W)])/
    (1+exp(param_vector[38+n.eta+2*(n.Z+n.W)]))
  
  a_sat      = bounds$min_a_sat+
    (bounds$max_a_sat-bounds$min_a_sat)*
    exp(param_vector[43+n.eta+2*(n.Z+n.W)])/
    (1+exp(param_vector[43+n.eta+2*(n.Z+n.W)]))
  
  sigma_H    = bounds$min_sigma_H+
    (bounds$max_sigma_H-bounds$min_sigma_H)*
    exp(param_vector[63+n.eta+2*(n.Z+n.W)])/
    (1+exp(param_vector[63+n.eta+2*(n.Z+n.W)]))
  
  new_model$parameters<-list(    
    A_bar      = A_bar,
    sigma_a    = sigma_a,
    Phi        = diag(exp(param_vector[3:(2+n.eta)])/
                        (1+exp(param_vector[3:(2+n.eta)]))),
    gamma      = exp(param_vector[2+n.eta+1]),
    delta      = exp(param_vector[2+n.eta+2]),
    m0         = exp(param_vector[2+n.eta+3]),
    ell0.N     = ell0.N,
    ell1.N     = ell1.N,
    rho.N      = rho.N,
    ell0.D     = ell0.D,                                                   
    ell1.D     = ell1.D,
    mu_d       = mu_d, 
    mu_n       = mu_n, 
    eps_0      = exp(param_vector[8+n.eta+2*(n.Z+n.W)+3]),                      #DICE2016,eland0
    rho        = exp(param_vector[8+n.eta+2*(n.Z+n.W)+4]),                      #DICE2016,deland
    mateq      = exp(param_vector[8+n.eta+2*(n.Z+n.W)+5]),                      #CDICE
    mueq       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+6]),                      #CDICE
    mleq       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+7]),                      #CDICE
    phi_0      = exp(param_vector[8+n.eta+2*(n.Z+n.W)+8]),                      #DICE2016,fex0
    phi_1      = exp(param_vector[8+n.eta+2*(n.Z+n.W)+9]),                      #DICE2016,fex1
    m_pi       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+10]),                     #CDICE
    xi_1       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+11]),                     #CDICE,c1
    xi_2       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+12]),                     #CDICE,c3
    xi_3       = exp(param_vector[8+n.eta+2*(n.Z+n.W)+13]),                     #CDICE,c4
    lambda     = exp(param_vector[8+n.eta+2*(n.Z+n.W)+14]),                     #DICE2016,noname=fco22x/t2xco2
    tau        = exp(param_vector[8+n.eta+2*(n.Z+n.W)+15]),                     #DICE2016,fco22x
    delsigma   =-exp(param_vector[8+n.eta+2*(n.Z+n.W)+16]),                     #DICE2016,dsig
    e0         = exp(param_vector[8+n.eta+2*(n.Z+n.W)+17]),                     #DICE2016
    q0         = exp(param_vector[8+n.eta+2*(n.Z+n.W)+18]),                     #DICE2016
    mu0        = exp(param_vector[8+n.eta+2*(n.Z+n.W)+19]),                     #DICE2016
    sigma0     = exp(param_vector[8+n.eta+2*(n.Z+n.W)+20]),                     #DICE2016
    gsigma1    =-exp(param_vector[8+n.eta+2*(n.Z+n.W)+21]),                     #DICE2016
    varphi_11  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+22]),                     #CDICE,b11
    varphi_12  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+23]),                     #CDICE,b12
    varphi_21  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+24]),                     #CDICE,b21
    varphi_22  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+25]),                     #CDICE,b22
    varphi_23  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+26]),                     #CDICE,b23
    varphi_32  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+27]),                     #CDICE,b32
    varphi_33  = exp(param_vector[8+n.eta+2*(n.Z+n.W)+28]),                     #CDICE,b33
    gback      = exp(param_vector[8+n.eta+2*(n.Z+n.W)+29]),                     #DICE2016
    pback      = pback,
    theta2     = exp(param_vector[8+n.eta+2*(n.Z+n.W)+31]),                     #DICE2016
    sigma_eta_f= exp(param_vector[8+n.eta+2*(n.Z+n.W)+32]),
    weights_Tg = exp(matrix(param_vector[(41+n.eta+2*(n.Z+n.W)):
                                           (42+n.eta+2*(n.Z+n.W))]
                            ,ncol=1)),
    a_sat      = a_sat,
    T_0s       = exp(param_vector[44+n.eta+2*(n.Z+n.W)]),
    b_sat      = exp(param_vector[45+n.eta+2*(n.Z+n.W)]),
    delta_K    = exp(param_vector[46+n.eta+2*(n.Z+n.W)]),
    ini_delc   = exp(param_vector[47+n.eta+2*(n.Z+n.W)]),
    ini_tildey = exp(param_vector[48+n.eta+2*(n.Z+n.W)]),
    ini_E      = exp(param_vector[49+n.eta+2*(n.Z+n.W)]),
    ini_Eind   = exp(param_vector[50+n.eta+2*(n.Z+n.W)]),
    ini_F      = exp(param_vector[51+n.eta+2*(n.Z+n.W)]),
    ini_Mat    = exp(param_vector[52+n.eta+2*(n.Z+n.W)]),
    ini_Mup    = exp(param_vector[53+n.eta+2*(n.Z+n.W)]),
    ini_Mlo    = exp(param_vector[54+n.eta+2*(n.Z+n.W)]),
    ini_Tat    = exp(param_vector[55+n.eta+2*(n.Z+n.W)]),
    ini_Tlo    = exp(param_vector[56+n.eta+2*(n.Z+n.W)]),
    ini_CumD   = exp(param_vector[57+n.eta+2*(n.Z+n.W)]),
    ini_CumE   = exp(param_vector[58+n.eta+2*(n.Z+n.W)]),
    ini_Cumdelc= exp(param_vector[59+n.eta+2*(n.Z+n.W)]),
    ini_H      = exp(param_vector[60+n.eta+2*(n.Z+n.W)]),
    tol.GN     = exp(param_vector[61+n.eta+2*(n.Z+n.W)]),
    eps.GN     = exp(param_vector[62+n.eta+2*(n.Z+n.W)]),
    sigma_H    = sigma_H,
    ell0.div   = exp(param_vector[64+n.eta+2*(n.Z+n.W)]),
    ell1.div   = matrix(exp(param_vector[(64+n.eta+2*(n.Z+n.W)+1):
                                    (64+n.eta+3*(n.Z+n.W))]),ncol=1),
    mu_div     = exp(param_vector[65+n.eta+3*(n.Z+n.W)]),
    chi.div    = exp(param_vector[66+n.eta+3*(n.Z+n.W)]),
    b_sk       = exp(param_vector[67+n.eta+3*(n.Z+n.W)]),
    c_sat      =-exp(param_vector[68+n.eta+3*(n.Z+n.W)])
  )
  
  new_model$ell1.N.tilde = ell1.N.tilde
  
  return(new_model)
}

#Loss function
#h is the period up to our optimization
compute.moments    <- function(param,
                         model,
                         FILTER,
                         horiz,
                         theta){
  h     <- horiz
  bounds<-model$bounds
  
  omega_ZCB     <-matrix(0,model_sol$n.X,1)
  param_vector<- Model2Param(model)
  if(length(param)!=sum(FILTER)){
    print("Problem length(param)!=sum(FILTER)")
    return(10^10)
  }
  param_vector[FILTER==1] <- param
  new_model <- Param2Model(param_vector,model)
  model_sol <- model_solve(new_model,theta)
  
  new_model_permafrost<- new_model
  model_sol_permafrost<- model_solve(new_model_permafrost,theta,
                                     indic_mitig=FALSE)

  new_model_nopermafrost <- new_model
  new_model_nopermafrost$parameters$mu_n <- 0
  model_sol_nopermafrost <- model_solve(new_model_nopermafrost,theta,
                                        indic_mitig=FALSE)
  
  EV              <- EV.fct(model_sol,h+40) ##!!
  EV_permafrost   <- EV.fct(model_sol_permafrost,h)
  EV_nopermafrost <- EV.fct(model_sol_nopermafrost,h)
  
  # Compute model-implied moments:
  model_implied_moment    <- NULL
  
  # Temperatures in 2100
  model_implied_moment[1] <- EV$EX[[9]][h]
  model_implied_moment[2] <- sqrt(EV$VX[[9]][h])
  
  # Diff in Temperatures wrt no permafrost:
  model_implied_moment[3] <- 0.2#EV_permafrost$EX[[9]][h]-EV_nopermafrost$EX[[9]][h]
  
  model_implied_moment[4] <- EV_permafrost$EX[[12]][h]-EV_nopermafrost$EX[[12]][h]
  
  # Slope Cum_D wrt Temp:
  #positive correlation (Cum_D, T_at)
  
  CovGMST<-EV_permafrost$CovX[[h]][11,9]
  VarGMST<-EV_permafrost$VX[[9]][h]
  
  model_implied_moment[5] <-CovGMST/VarGMST
  
  # Very-long-term rate
  v.0                     <-varphi(model_sol,model_sol$omega_ZCB,model_sol$Tmax)                                              
  model_implied_moment[6] <-v.0$r.t[h]                                          #LT interest rate, in %
  
  # sea level targets
  model_implied_moment[7] <- EV$EX[[14]][h]
  model_implied_moment[8] <- sqrt(EV$VX[[14]][h])
  
  # non explosive model
  model_implied_moment[9]  <- EV$EX[[5]][h]
  model_implied_moment[10] <- EV$EX[[5]][h+20]
  model_implied_moment[11] <- EV$EX[[5]][h+40]
  
  model_implied_moment[12] <- EV$EX[[9]][h+20]
  model_implied_moment[13] <- EV$EX[[9]][h+40]
  
  model_implied_moment[14] <- min(v.0$r.t[1:model_sol$Tmax]) 
  
  # Temperature in 2040-2050
  model_implied_moment[15] <- EV$EX[[9]][h-12]
  model_implied_moment[16] <- EV$EX[[9]][h-10]
  
  
  # return crazy value when NaNs have been obtained in computation of mu.u0.t
  if(model_sol$u0 == 10000){
    model_implied_moment <- rep(100,length(model_implied_moment))
  }
  return(model_implied_moment)
}

lossfunction4optim <- function(param,
                         model,
                         FILTER,
                         horiz,
                         theta){
  model_implied_moment<-compute.moments(param,
                                        model,
                                        FILTER,
                                        horiz,
                                        theta)
  weights      <-model$weights
  #print(model_implied_moment)
  
  # Compute loss function:
  loss <- sum(weights * (model_implied_moment - target_vector)^2
              +(model_implied_moment[13]>model_implied_moment[12])*
                (model_implied_moment[13]-model_implied_moment[12])^2*100
              +(model_implied_moment[14]<0)*
                (model_implied_moment[14]-0)^2*1000
  )
  return(loss)
}


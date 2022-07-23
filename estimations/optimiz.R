#Functions
source("./procedure/Functions_General_v2.R")
source("./procedure/Functions_Optim.R")
# =========================================================
# Step 1 - Estimate parameters except Sea Level ones:
# Maximum number of iterations for Nelder-Mead:
MAXIT.NldMd <- 150
# Maximum number of iterations for nlminb:
MAXIT.nlminb<- 0
# Number of loops:
nb.loops    <- 4

#* FILTER is a vector of binary entries (0 or 1). If i^th component is 1, then
#*  this component can be changed by optimization:
FILTER     <- rep(0,68+n.eta+3*(n.Z+n.W))                                       #CHANGE nb param.
FILTER[1]                     <- 1 # A_bar
FILTER[2]                     <- 1 # sigma_a
FILTER[6+n.eta]               <- 1 # ell.0.N
FILTER[15+n.eta]              <- 1 # ell.1.N[9]
FILTER[7+n.eta+n.Z+n.W]       <- 1 # rho.N
FILTER[17+n.eta+n.Z+n.W]      <- 1 # ell.1.D[9]
FILTER[9+n.eta+2*(n.Z+n.W)]   <- 1 # mu.D
FILTER[10+n.eta+2*(n.Z+n.W)]  <- 1 # mu.N
FILTER[38+n.eta+2*(n.Z+n.W)]  <- 1 # pback
FILTER[40+n.eta+2*(n.Z+n.W)]  <- 1 # sigma_eta_f


# Initialize estimated model:
model_est <- model

if(indic_estim==1){
  print("Starting optimization: takes around 30 minutes.")
  print("4 loops of Nelder-Mead optimization... I might take few minutes.")
  source("estimations/aux_use_Numerical_algo.R")
}else{
  load(file="estimations/solve_param_estim/param_final_est.Rdat")
}

# =========================================================
# Step 2 - Estimate Sea Level parameters:
# Maximum number of iterations for Nelder-Mead:
MAXIT.NldMd <- 150
# Maximum number of iterations for nlminb:
MAXIT.nlminb<- 0
# Number of loops:
nb.loops    <- 2

#RESET FILTER
FILTER                        <- rep(0,length(FILTER))
FILTER[43+n.eta+2*(n.Z+n.W)]  <- 1 # a_sat
FILTER[68+n.eta+3*(n.Z+n.W)]  <- 1 # c_sat

if(indic_estim==1){
  print("2 loops of Nelder-Mead optimization (sea level)... I might take few minutes.")
  source("estimations/aux_use_Numerical_algo.R")
}else{
  load(file="estimations/solve_param_estim/param_final_est.Rdat")
}

# =========================================================
# Step 3 - Test:
# Maximum number of iterations for Nelder-Mead:
MAXIT.NldMd <- 150
# Maximum number of iterations for nlminb:
MAXIT.nlminb<- 10
# Number of loops:
nb.loops    <- 2

#RESET FILTER
FILTER                        <- rep(0,length(FILTER))
FILTER[1]                     <- 1 # A_bar
FILTER[2]                     <- 1 # sigma_a
FILTER[6+n.eta]               <- 1 # ell.0.N
FILTER[15+n.eta]              <- 1 # ell.1.N[9]
FILTER[7+n.eta+n.Z+n.W]       <- 1 # rho.N
FILTER[17+n.eta+n.Z+n.W]      <- 1 # ell.1.D[9]
FILTER[9+n.eta+2*(n.Z+n.W)]   <- 1 # mu.D
FILTER[10+n.eta+2*(n.Z+n.W)]  <- 1 # mu.N
FILTER[38+n.eta+2*(n.Z+n.W)]  <- 1 # pback
FILTER[40+n.eta+2*(n.Z+n.W)]  <- 1 # sigma_eta_f
FILTER[43+n.eta+2*(n.Z+n.W)]  <- 1 # a_sat
FILTER[68+n.eta+3*(n.Z+n.W)]  <- 1 # c_sat

if(indic_estim==1){
  print("2 loops of Nelder-Mead + Nlminb optimizations (all parameters)... I might take few minutes.")
  source("estimations/aux_use_Numerical_algo.R")
}else{
  load(file="estimations/solve_param_estim/param_final_est.Rdat")
}


# Retrieve estimated parameters:
param_vector <- Model2Param(model_est)
param_est    <- param_vector[FILTER==1]

#Show only 4 decimals
nb.dec   <-4
format.nb<-paste("%.",nb.dec,"f",sep="")
#sprintf(format.nb,1)


model_sol<-model_solve(model_est,theta0)

mim<-compute.moments(param_est,model_est,FILTER,horiz,theta0)
print("New moments and targets achieved")
print(cbind(target_vector,mim))
print(cbind(c("mu_n","mu_d","rho.N","ell0.N","ell1.N","ell1.D","pback",
              "a_sat","c_sat"),
            c(sprintf(format.nb,model_est$parameters$mu_n),
              sprintf(format.nb,model_est$parameters$mu_d),
              sprintf(format.nb,model_est$parameters$rho.N),
              sprintf(format.nb,model_est$parameters$ell0.N),
              sprintf(format.nb,model_est$parameters$ell1.N[9]),
              sprintf(format.nb,model_est$parameters$ell1.D[9]),
              sprintf(format.nb,model_est$parameters$pback),
              sprintf(format.nb,model_est$parameters$a_sat),
              sprintf(format.nb,model_est$parameters$c_sat)
              )))
tic("Loss Function")
min_loss<-lossfunction4optim(param_est,
                   model_est,
                   FILTER,
                   horiz=horiz,
                   theta0)
toc()
print("Minimized Loss")
print(min_loss)

print("Loss decomposition")
print(weights * (mim - target_vector))
plot(mim[-4],col="red",ylim=c(range(mim[-4],target_vector[-4])))
lines(target_vector[-4],type="p")
abline(h=0,col="green")

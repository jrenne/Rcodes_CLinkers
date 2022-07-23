# Build vector of parameters corresponding to initial model:
param_vector_ini <- Model2Param(model_est)

# Extract the components of param_vector_ini that will be optimized on:
param_ini <- param_vector_ini[FILTER==1]
tic("Loss Function")
loss<-lossfunction4optim(param_ini,
                         model_est,
                         FILTER,
                         horiz=horiz,
                         theta0)
toc()

print("Initial Loss")
print(loss)

# param_est will be updated after each use of numerical optim algorithm:
param_est <- param_ini 
# number of passes in the two algorithms
all.criteria<-NULL


for(iiii in 1:nb.loops){
  RES.optim <- optimx(param_est,
                      lossfunction4optim,
                      model = model_est,
                      FILTER = FILTER,
                      horiz=horiz,
                      theta=theta0,
                      method = "Nelder-Mead",
                      control=list(trace=1, maxit = MAXIT.NldMd,
                                   kkt = FALSE))
  
  # To make sure that model is updated only if optim results are not NaN:
  if(sum(is.na(as.matrix(RES.optim)[1:length(param_est)]))==0){
    param_est     <- c(as.matrix(RES.optim)[1:length(param_est)])
    save(param_est,file="estimations/solve_param_estim/param_NldMd.Rdat")
    all.criteria<-c(all.criteria,RES.optim$value)
  }
  
  if(MAXIT.nlminb>0){
    RES.optim <- optimx(param_est,
                        lossfunction4optim,
                        model = model_est,
                        FILTER = FILTER,
                        horiz=horiz,
                        theta=theta0,
                        method = "nlminb",
                        control=list(trace=1, maxit = MAXIT.nlminb,
                                     kkt = FALSE))
    
    # To make sure that model is updated only if optim results are not NaN:
    if(sum(is.na(as.matrix(RES.optim)[1:length(param_est)]))==0){
      param_est     <- c(as.matrix(RES.optim)[1:length(param_est)])
      save(param_est,file="estimations/solve_param_estim/param_nlminb.Rdat")
      all.criteria<-c(all.criteria,RES.optim$value)
    }
  }
}

# Save estimated model:
param_vector_est            <- Model2Param(model_est)
param_vector_est[FILTER==1] <- param_est
model_est                   <- Param2Model(param_vector_est,model_est)
save(model_est,file="estimations/solve_param_estim/param_final_est.Rdat")


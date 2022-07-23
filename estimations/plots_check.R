#Verification plots
source("./procedure/Functions_General_v2.R")

#Simulations
nb.simul<-model_sol$Tmax-1                                                                   
nb.traj <-10000                                                                 
H       <-model$horiz.2100



x<-exp(seq(-20,40,length.out = 1000))                                           #grid for proposition 8 (Fourier)
a<-omega_T.at                                                                   #Options related to T_at in X
b<-2                                                                            #Options pay when a'X < b

##################
#Model-implied EV#
##################
EV<-EV.fct(model_sol)

##########################
#Compare with simulations#
##########################
test <-simul.function(model_sol,nb.simul,nb.traj)
test <-test[-length(test)] #vectors of X

##Create plots with mean and confidence intervals
#Quantiles
quantiles1<-c(0.025,0.975)
Q<-lapply(test,function(m)
  apply(m,1,function(x)quantile(x,quantiles1,na.rm=TRUE)))
#Variance
sample.var<-lapply(test,function(m)apply(m,1,var))
#Mean of simulations
mean<-lapply(test, function(m)apply(m,1,mean))[1:15]

sample.upper<-lapply(1:length(mean),function(x)mean[[x]]+2*sample.var[[x]]^(1/2))
sample.lower<-lapply(1:length(mean),function(x)mean[[x]]-2*sample.var[[x]]^(1/2))

bounds.simul<-lapply(1:length(mean),function(x){
  rbind(sample.lower[[x]],sample.upper[[x]])
})
remove(sample.upper, sample.lower)


#########
#mu plot#
#########
mudice<-0
mudice[1:(length(model_sol$vec_date)+1)]<-pmin(exp(-1/21*log(0.17)*
                                          (1:(length(model_sol$vec_date)+1)))*
                                            0.17,1)                             #mu from DICE
plot(model_sol$mu,type="l")
lines(mudice,col="red")
legend("bottomright",legend=c("DICE","Own mitigation"),col=c("red","black"),
       lty=c(1,1))


#########
#Pricing#
#########

T.swaps<-varphi.tilde(model_sol,omega_T.at,H)[[1]]/
   varphi(model_sol,omega_ZCB,H)[[3]]

print(T.swaps-EV$EX$T_at[1:H])


#############
#Optim Plots#
#############

#Check for target [5]

plot(test$T_at[H,],test$Cum_dc[H,],
     main="Correlation damages and temperatures")
abline(lm(test$Cum_dc[H,]~test$T_at[H,]))

#Check for targets [9:13]
non_explosive<-EV.fct(model_sol,200)
plot(non_explosive$date,non_explosive$EX$T_at,type="l")
mod<-model_solve(model_sol,theta0)

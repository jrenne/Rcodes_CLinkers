par(plt=c(.1,.9,.1,.9))

values.of.mu_d <- seq(0.6*model_sol$parameters$mu_d,
                      1.4*model_sol$parameters$mu_d,length.out = 8)
values.of.T_at <- seq(1,4,length.out=8)
# if(plots_run_all==1){
# }
all.models <- list()
all.EV     <- list()

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/parallel.Rdata")
clusterEvalQ(cl,load("outputs/parallel.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.Ts <- foreach(i = 1:length(values.of.mu_d), .combine=rbind) %dopar% {
  
  model_i <- model_sol
  
  model_i$parameters$mu_d <- values.of.mu_d[i]
  model_i_sol<-model_solve(model_i,theta0)
  all.models[[i]] <- model_i_sol
  
  model_i_P <- model_i_sol
  model_i_P$parameters$gamma <- 1.0000000001
  model_i_sol_P <- model_solve(model_i_P,theta0)
  
  EV.i       <- EV.fct(model_i_sol,h=H)
  all.EV[[i]]<- EV
  T.P        <- EV.i$EX$T_at[H]
  varphi.i   <- varphi.tilde(model_i_sol,omega_T.at,H)[[1]]/
    varphi(model_i_sol,omega_ZCB,H)[[3]]
  T.Q        <- varphi.i[H]
  # Check:
  varphi.i.P<- varphi.tilde(model_i_sol_P,omega_T.at,H)[[1]]/
    varphi(model_i_sol_P,omega_ZCB,H)[[3]]
  T.P.check <- varphi.i.P[H]
  
  c(T.P,T.Q,T.P.check)
}

stopCluster(cl)


T.P       <- matrix(all.Ts[,1],length(values.of.mu_d),1)
T.Q       <- matrix(all.Ts[,2],length(values.of.mu_d),1)
T.P.check <- matrix(all.Ts[,3],length(values.of.mu_d),1)

#Plots
#1
FILE = paste("/outputs/Figures/Figure_cut_CP.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=6)

par(mfrow=c(1,1))
par(plt=c(.15,.95,.15,.85))

# Refine grid (interpolation):
values.of.mu_d.fine <- seq(values.of.mu_d[1],
                           tail(values.of.mu_d,1),length.out = 40)

# Use splines:
spline.T.P<-splinefun(values.of.mu_d,T.P,method="natural")
spline.T.Q<-splinefun(values.of.mu_d,T.Q,method="natural")

T.P.fit <- spline.T.P(values.of.mu_d.fine)
T.Q.fit <- spline.T.Q(values.of.mu_d.fine)

plot(values.of.mu_d.fine,T.P.fit,type="l",
     col=P.col.line,lwd=2,las=1,
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     main="",
     xlab=expression(paste("Disaster magnitude ",mu[D],sep="")),
     ylab=expression(paste("Temperature Anomaly ",T[at],sep="")),
     ylim=c(range(T.P.fit,T.Q.fit)))
lines(values.of.mu_d.fine,T.Q.fit,
      col=Q.col.line,lwd=2)
abline(v=model_sol$parameters$mu_d,col="grey",lty=3,lwd=2)
legend("bottomleft",
       legend=c(expression(paste("Expected ",T[at]," in 2100",sep="")),
                expression(paste("Swaps Price ",T^S," in 2100",sep=""))
                ),
       lty=c(1,1),
       lwd=rep(2,2),
       col=c(P.col.line,Q.col.line),
       bty = "n",cex=1.2)
dev.off()


#Climate beta
y.lim   <- c(0,2)
y.limscc<- c(0,2000)
y.limzcb<- c(-4,5)
cexaxs  <- 2
cexlab  <- 2.5
cexmain <- 2.5

# Maturity, in periods, of the ZCB:
zcbmat  <-10

# Considered values of mu_n (select three values):
values.of.mu_n <- c(0,
                    model_sol$parameters$mu_n,
                    40)

# Considered values of mu_d:
seq.mu_d<-seq(0,0.15,len=80)


# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("outputs/parallel.Rdata")
clusterEvalQ(cl,load("outputs/parallel.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
print("Climate Beta figure -- Loop on mu_d values")
all_res <- foreach(i = 1:length(seq.mu_d), 
                   .combine=rbind) %dopar% {
                     
                     modelswp1<-model_sol
                     modelswp2<-model_sol
                     modelswp3<-model_sol
                     
                     modelswp1$parameters$mu_d<-seq.mu_d[i]
                     modelswp2$parameters$mu_d<-seq.mu_d[i]
                     modelswp3$parameters$mu_d<-seq.mu_d[i]
                     
                     modelswp1$parameters$mu_n <- values.of.mu_n[1]
                     modelswp2$parameters$mu_n <- values.of.mu_n[2]
                     modelswp3$parameters$mu_n <- values.of.mu_n[3]
                     
                     model_solswp1<-model_solve(modelswp1,theta0)
                     model_solswp2<-model_solve(modelswp2,theta0)
                     model_solswp3<-model_solve(modelswp3,theta0)
                     
                     # Get values for model1:
                     swaps1    <-(varphi.tilde(model_solswp1,omega_T.at,H)[[1]]/
                                    varphi(model_solswp1,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp1,H)$EX$T_at[H]
                     scc1      <- scc.fct(model_solswp1,0)
                     zcb1      <-varphi(model_solswp1,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     # Get values for model2:
                     swaps2    <-(varphi.tilde(model_solswp2,omega_T.at,H)[[1]]/
                                    varphi(model_solswp2,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp2,H)$EX$T_at[H]
                     scc2      <- scc.fct(model_solswp2,0)
                     zcb2      <-varphi(model_solswp2,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     # Get values for model3:
                     swaps3    <-(varphi.tilde(model_solswp3,omega_T.at,H)[[1]]/
                                    varphi(model_solswp3,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp3,H)$EX$T_at[H]
                     scc3      <- scc.fct(model_solswp3,0)
                     zcb3      <-varphi(model_solswp3,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     
                     c(swaps1,swaps2,swaps3,scc1,scc2,scc3,zcb1,zcb2,zcb3)
                   }
print("End of computations for Climate Beta figure")

stopCluster(cl)

all_res[abs(all_res)>10^5] <- NaN

swaps     <-all_res[,2]
scc       <-all_res[,5]
zcb       <-all_res[,8]

swaps0    <-all_res[,1]
scc0      <-all_res[,4]
zcb0      <-all_res[,7]

swaps50   <-all_res[,3]
scc50     <-all_res[,6]
zcb50     <-all_res[,9]

splineswaps50 <-splinefun(seq.mu_d[!is.na(swaps50)],swaps50[!is.na(swaps50)],method = "natural")
splineswaps0  <-splinefun(seq.mu_d[!is.na(swaps0)],swaps0[!is.na(swaps0)],method = "natural")
splineswaps   <-splinefun(seq.mu_d[!is.na(swaps)],swaps[!is.na(swaps)],method = "natural")

splinescc50   <-splinefun(seq.mu_d[!is.na(scc50)],scc50[!is.na(scc50)],method = "natural")
splinescc0    <-splinefun(seq.mu_d[!is.na(scc0)],scc0[!is.na(scc0)],method = "natural")
splinescc     <-splinefun(seq.mu_d[!is.na(scc)],scc[!is.na(scc)],method = "natural")

splinezcb50   <-splinefun(seq.mu_d[!is.na(zcb50)],zcb50[!is.na(zcb50)],method = "natural")
splinezcb0    <-splinefun(seq.mu_d[!is.na(zcb0)],zcb0[!is.na(zcb0)],method = "natural")
splinezcb     <-splinefun(seq.mu_d[!is.na(zcb)],zcb[!is.na(zcb)],method = "natural")

grid          <-seq(seq.mu_d[1],tail(seq.mu_d,1),length.out = 200+1)

swaps.50<-splineswaps50(grid)
swaps.0 <-splineswaps0(grid)
swaps.n <-splineswaps(grid)

scc.50  <-splinescc50(grid)
scc.0   <-splinescc0(grid)
scc.n   <-splinescc(grid)

zcb.50  <-splinezcb50(grid)
zcb.0   <-splinezcb0(grid)
zcb.n   <-splinezcb(grid)

climate.beta<-cbind(grid,swaps.n,swaps.50,swaps.0,
                    scc.n,scc.50,scc.0,
                    zcb.n,zcb.50,zcb.0)
save(climate.beta,file="./outputs/backup_parallel/climate_data.Rdat")

#Plots
mun<-as.integer(model_sol$parameters$mu_n)

FILE = paste("/outputs/Figures/Figure_Climate_Beta.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=10, height=10)
par(plt=c(.1,.95,.1,.9))
par(mfrow=c(3,1))

plot(grid,swaps.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.lim,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(mu[D]), main="(a) - Risk-Adjusted Temperature in 2100",
     ylab=expression(paste("Swap Price minus Expected ", T[AT*','*2100]),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,swaps.0,col="tan1",
      lty=2,lwd=2)
lines(grid,swaps.50,col="firebrick1",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (calibration)"),
                          expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=2.5)

plot(grid,scc.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.limscc,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(mu[D]),main="(b) - Social Cost of Carbon",
     ylab=expression(paste("SCC (in $/GTC)"),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,scc.0,col="tan1",
      lty=2,lwd=2)
lines(grid,scc.50,col="firebrick1",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (calibration)"),
                          expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=2.5) 
plot(grid,zcb.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.limzcb,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(mu[D]),main="(c) - Long-term Real Interest Rate (50 years)",
     ylab=expression(paste("Real Rate (in percent)"),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,zcb.0,col="tan1",
      lty=2,lwd=2)
lines(grid,zcb.50,col="firebrick1",
      lty=4,lwd=2)
legend("bottomleft",legend=c(expression(paste(mu[N],"=0")),
                             bquote(mu[N]==.(mun) ~" (calibration)"),
                             expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=2.5)
dev.off()

#Separated Plots
###
FILE = paste("/outputs/Figures/Figure_Climate_Beta_swaps.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=6)
par(mfrow=c(1,1))
plot(grid,swaps.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.lim,
     xlab=expression(mu[D]), main="Risk-Adjusted Temperature in 2100",
     ylab=expression(paste("Swap Price minus Expected ", T[AT*','*2100]),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,swaps.0,col="tan1",
      lty=2,lwd=2)
lines(grid,swaps.50,col="firebrick1",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (calibration)"),
                          expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()
##
FILE = paste("/outputs/Figures/Figure_Climate_Beta_SCC.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=6)
par(mfrow=c(1,1))
plot(grid,scc.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.limscc,
     xlab=expression(mu[D]),main="Social Cost of Carbon",
     ylab=expression(paste("SCC (in $/GTC)"),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,scc.0,col="tan1",
      lty=2,lwd=2)
lines(grid,scc.50,col="firebrick1",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (calibration)"),
                          expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()
#
FILE = paste("/outputs/Figures/Figure_Climate_Beta_ZCB.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=6)
par(mfrow=c(1,1))
plot(grid,zcb.n,type="l",lwd=2,col=Q.col.line,
     ylim=y.limzcb,
     xlab=expression(mu[D]),main="Long-term Real Interest Rate (50 years)",
     ylab=expression(paste("Real Rate (in percent)"),
                     sep=""),
     las=1)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_d,lty=3,col="lightslategrey")
lines(grid,zcb.0,col="tan1",
      lty=2,lwd=2)
lines(grid,zcb.50,col="firebrick1",
      lty=4,lwd=2)
legend("bottomleft",legend=c(expression(paste(mu[N],"=0")),
                             bquote(mu[N]==.(mun) ~" (calibration)"),
                             expression(paste(mu[N],"=40"))),
       col=c("tan1",Q.col.line,"firebrick1"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()
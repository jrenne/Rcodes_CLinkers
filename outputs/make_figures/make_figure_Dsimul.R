#Gamma 0 shock simulations
nb      <-12
nb.simul<-model_sol$Tmax                                                                   

simul_d <-simul.function(model_sol,nb.simul,nb)

D <-simul_d$Cum_D[1:H,]
N <-simul_d$N[1:H,]

Cum_N<-matrix(1,H,nb)

while(min(Cum_N[H,])!=0){
  simul_d <-simul.function(model_sol,nb.simul,nb)
  
  D <-simul_d$Cum_D[1:H,]
  N <-simul_d$N[1:H,]
  
  Cum_N<-matrix(N[1,],H,nb)
  for(i in 1:(H-1)){
    Cum_N[i+1,]<-Cum_N[i,]+N[i+1,]
  } 
}

ED<-EV$EX$Cum_D
EN<-EV$EX$N
Cum_EN<-EV$EX$N[1]
for(i in 1:(H-1)){
  Cum_EN[i+1]<-Cum_EN[i]+EN[i+1]
}

ED<-exp(ED)-1
D <-exp(D) -1

#bigN  <-Cum_N[,max.col(Cum_N)[H]]
smallN<-Cum_N[,max.col(-Cum_N)[H]]
#bigD  <-D[,max.col(-D)[H]]
smallD<-D[,max.col(D)[H]]

#Plots
#D
FILE = "/outputs/Figures/Figure_Cum_Damages_simul.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=10, height=8)
par(mfrow=c(2,1))  
plot(vec_date[2:(H+1)],ED,col=P.col.line,type="l",
     ylim=c(range(ED,D)),
     lwd=2,lty=2,las=1,
     ylab="Consumption Damages (in percent)", xlab="",
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
for(i in 1:dim(D)[2]){
  lines(vec_date[2:(H+1)],D[,i],
        col="grey75")
}
lines(vec_date[2:(H+1)],smallD,
      col=P.col.line)
legend("bottomleft",legend=c("Expected Path","Simulations"),
       col=c(P.col.line,"grey75"),lty=c(2,1),lwd=c(2,1),
       bty = "n",cex=1.2)
#N
plot(vec_date[2:(H+1)],Cum_EN,col=P.col.line,type="l",
     ylim=c(range(Cum_EN,Cum_N)),
     lwd=2,lty=2,las=1,
     ylab="Emissions from Permafrost", xlab="",
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
for(i in 1:dim(Cum_N)[2]){
  lines(vec_date[2:(H+1)],Cum_N[,i],
        col="grey75")
}
lines(vec_date[2:(H+1)],smallN,
      col=P.col.line)
legend("topleft",legend=c("Expected Path","Simulations"),
       col=c(P.col.line,"grey75"),lty=c(2,1),lwd=c(2,1),
       bty = "n",cex=1.2)

dev.off()

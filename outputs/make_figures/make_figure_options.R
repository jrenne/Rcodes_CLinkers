#Options on temperatures (atm.)
X.tat<-9

q3 <- sequential_hcl(5, "YlOrRd")
q3 <-q3[3:1]

x.lim  <-c(2020,2100)
pos.opt<-c(3,3,3)

P.OpD.names<-c(expression(paste("T"[K],"=2.5"),
                          paste("T"[K],"=3"),
                          paste("T"[K],"=4")))


P.OpD     <-list()
P.OpD_RNF <-list(1,1,1)                                                         #Risk-neutrality w/ Fourier

#Normalization with ZCB
v.0<-varphi(model_sol,omega_ZCB,H) 


for(k in 1:length(K)){
  v.hat     <-varphi.hat(model_sol,omega_ZCB,H,x,-a,-K[k])
  P.hat.N   <-v.hat$P.hat/v.0$P.t
  P.OpD[[k]]<-P.hat.N*100
}
#Fourier
for(k in 1:length(K)){
  for( i in 1:H){
    P.OpD_RNF[[k]][i]<-(1-fourier(model_sol,x,K[k],i,X.tat))*100
  }
}


for(k in 1:length(K)){
  P.OpD_RNF[[k]][1]<-P.OpD_RNF[[k]][2]
  P.OpD[[k]][1]    <-P.OpD[[k]][2]
  for(i in 3:H){
    if(P.OpD_RNF[[k]][i]< P.OpD_RNF[[k]][i-1]){
      P.OpD_RNF[[k]][i] <-P.OpD_RNF[[k]][i-1]
    }
    if(P.OpD[[k]][i]< P.OpD[[k]][i-1]){
      P.OpD[[k]][i] <-P.OpD[[k]][i-1]
    }
  }
}


splinefunP1    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD_RNF[[1]],
                           method="hyman")
splinefunQ1    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD[[1]],
                           method="hyman")
splinefunP2    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD_RNF[[2]],
                           method="hyman")
splinefunQ2    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD[[2]],
                           method="hyman")
splinefunP3    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD_RNF[[3]],
                           method="hyman")
splinefunQ3    <-splinefun(model_sol$vec_date[2:(H+1)],P.OpD[[3]],
                           method="hyman")
P.OpDPs<-list()
P.OpDQs<-list()
grid   <-model_sol$vec_date[2]:model_sol$vec_date[H+1]

P.OpDPs[[1]]<-splinefunP1(grid)
P.OpDPs[[2]]<-splinefunP2(grid)
P.OpDPs[[3]]<-splinefunP3(grid)

P.OpDQs[[1]]<-splinefunQ1(grid)
P.OpDQs[[2]]<-splinefunQ2(grid)
P.OpDQs[[3]]<-splinefunQ3(grid)

text.pos  <-c(model_sol$vec_date[2]+20-2,
              model_sol$vec_date[2]+40-2,
              model_sol$vec_date[2]+60-2)
text.pos.y<-c(P.OpDQs[[1]][grep(text.pos[1]+2,model_sol$vec_date)]+1,
              P.OpDQs[[2]][grep(text.pos[2]+2,model_sol$vec_date)]+1,
              P.OpDQs[[3]][grep(text.pos[3]+2,model_sol$vec_date)]+1)

#Plot
FILE = "/outputs/Figures/Figure_Option_Digital.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=10, height=6)
par(mar=c(5.1,4.1,3,2.1))
plot(grid,P.OpDPs[[1]],type="l",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",
     lwd=2,lty=2,
     ylab="In percent",xlab="Maturity",
     las=1,
     ylim=c(range(P.OpD,P.OpD_RNF)),xlim=x.lim)
for (k in 1:length(K)){
  lines(grid,P.OpDQs[[k]],
        col=q3[k],
        lwd=2,
        lty=1)
  lines(grid,P.OpDPs[[k]],
        col=q3[k],
        lwd=2,
        lty=2)
  text(text.pos[k],text.pos.y[k],
       labels = P.OpD.names[k],
       cex = 1.2, pos = pos.opt[k], col = q3[k])
}
legend("topleft",
       legend=c("Option price","(= risk-adjusted probability)",
                "Option price without risk premium","(= physical probability)"),
       lty=c(1,0,2),
       lwd=rep(2,0,2),
       col=c("black","white","black","white"),
       bty = "n",cex=1.2)
dev.off()

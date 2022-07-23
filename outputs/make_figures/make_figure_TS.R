#Swaps and TIBs

tib<-list()
chi<-c(0.5,1,2)
for(i in 1:length(chi)){
  tib[[i]]<-TIB(model_sol,chi[i],EV$EX[[9]][1:H],H,9)
}
lty.tib<-c(2,1,4)
chi.names<-c(expression(paste(chi,"=0.5")),
             expression(paste(chi,"=1")),
             expression(paste(chi,"=2")))

T.swaps<-varphi.tilde(model_sol,omega_T.at,H)[[1]]/
          varphi(model_sol,omega_ZCB,H)[[3]]
v.0    <-varphi(model_sol,omega_ZCB,H)

FILE = "/outputs/Figures/Figure_TIB_Swaps.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=10, height=8)
par(mfrow=c(2,1))
plot(vec_date[2:(H+1)],EV$EX$T_at[1:H],
     type="l",las=1,
     xlab="Maturity",
     ylab="Temperature Anomaly",
     main="(a) - Term structures of Temperature-Indexed Swap rates",
     ylim = c(range(EV$EX$T_at[1:H],T.swaps)),
     lwd=2,
     lty=1,
     col=P.col.line,
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
#cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(vec_date[2:(H+1)],T.swaps,
      col=Q.col.line,
      lwd=2)
legend("topleft",
       legend=c("Expected Temperature Path",
                expression(paste("Swaps Prices T"^"S"))),
       lty=c(1,1),
       col=c(P.col.line,Q.col.line),
       lwd=c(2,2),
       bty = "n",cex=1.2)
#cex=1.5)
#par(mar=c(5.1,4.1,4.1,2.1))
plot(vec_date[2:(H+1)],v.0$r.t,
     type="l",
     ylim = c(range(tib[[3]][[2]],v.0$r.t)),
     #xlim=x.lim,
     main="(b) - Term structures of Temperature-Indexed Bonds yields",
     col="black",
     lwd=2,
     xlab = "Maturity",
     ylab="Interest Rates (in percent)",
     las=1,
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
for(i in 1:length(chi)){
  lines(vec_date[2:(H+1)],tib[[i]]$r.tib,
        col=Q.col.line,
        lty=lty.tib[i],
        lwd=2)
}
legend("bottomleft",legend=c("TIB","ZCB"),col=c(Q.col.line,"black"),
       lty =c(1,1),
       lwd=c(2,2),
       bty = "n",cex=1.2)
legend("topright",legend=c(chi.names[1],chi.names[2],chi.names[3]),
       col=c(Q.col.line,Q.col.line,Q.col.line),
       lty =c(2,1,4),
       lwd=c(2,2,2),
       text.col=c(Q.col.line,Q.col.line,Q.col.line),
       bty = "n",cex=1.2)
dev.off()


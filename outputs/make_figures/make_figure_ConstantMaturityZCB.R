indic_stocks<-0
# 50-Y constant maturity
h_cst<-10
h_cst.2<-2
# over a 50-Y window
nb<-model_sol$horiz.2100+1

# Model with climate change (CC) - mu_d >0, b_sk >0
model_cc <-model_sol

cc.10<-compute_cst_h(model_cc,h_cst,nb)
cc.2 <-compute_cst_h(model_cc,h_cst.2,nb)


# Model without climate change risk but agent does not re-optimize
model_ncc<-model_cc
model_ncc$parameters$mu_d<-0
model_ncc$parameters$b_sk<-0

model_ncc<-model_solve(model_ncc,theta0,indic_mitig = F,mu.chosen=model_sol$mu)

ncc.10<-compute_cst_h(model_ncc,h_cst,nb)
ncc.2 <-compute_cst_h(model_ncc,h_cst.2,nb)

 
# Model without climate change risk + agent re-optimizes
model_ncc.opt<-model_solve(model_ncc,theta0)

ncc.opt.10<-compute_cst_h(model_ncc.opt,h_cst,nb)
ncc.opt.2 <-compute_cst_h(model_ncc.opt,h_cst.2,nb)


# Model without climate change risk + agent re-optimizes + mu_c=E_t(delc)
indic_tp<-1

model_ncc.mu_c     <- model_ncc
model_ncc.mu_c$mu_c<- c(model_cc$X[1],EV.fct(model_cc)$EX$delc)
model_ncc.mu_c     <- model_solve(model_ncc.mu_c,theta0)

ncc.mu_c.10<-compute_cst_h(model_ncc.mu_c,h_cst,nb)
ncc.mu_c.2 <-compute_cst_h(model_ncc.mu_c,h_cst.2,nb)

indic_tp<-0

# 4 models
col<-c(brocolors("crayons")["Orange Red"],
       brocolors("crayons")["Wild Blue Yonder"],
       brocolors("crayons")["Indigo"],
       brocolors("crayons")["Violet Blue"])


#Plot
cexaxs  <- 2
cexlab  <- 2.5
cexmain <- 2.5

FILE = paste("/outputs/Figures/Figure_cstH_ZCB.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=9, height=9)
par(plt=c(.1,.95,.1,.9))
par(mfrow=c(2,1))
plot(cc.10[,1],cc.10[,2],
     type="l", xlab="Date",ylab="Interest rate (in percent)",
     ylim=c(range(cc.10[,2]-0.5,ncc.10[,2],ncc.opt.10[,2]+0.2,ncc.mu_c.10[,2])),
     las=1,
     cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
     main="50-Y constant maturity interest rate",
     lwd=2,col=col[1])
lines(ncc.10[,1],ncc.10[,2],
      lwd=2,col=col[2])
lines(ncc.opt.10[,1],ncc.opt.10[,2],
      lwd=2,col=col[3])
lines(ncc.mu_c.10[,1],ncc.mu_c.10[,2],
      lwd=2,col=col[4])
legend("bottomleft",
       legend=c("w/ CC (Climate Change)",
                "w/o CC + w/o optim.",
                "w/o CC + w/ optim.",
                "w/o CC + w/ optim. + imposed growth path"),
       lty=rep(1,4),
       col=col,
       lwd=rep(2,4),
       bty = "n",cex=2)

plot(cc.2[,1],cc.2[,2],
     type="l", xlab="Date",ylab="Interest rate (in percent)",
     ylim=c(range(cc.2[,2]-0.5,ncc.2[,2],ncc.opt.2[,2]+0.2,ncc.mu_c.2[,2])),
     las=1,
     cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
     main="10-Y constant maturity interest rate",
     lwd=2,col=col[1])
lines(ncc.2[,1],ncc.2[,2],
      lwd=2,col=col[2])
lines(ncc.opt.2[,1],ncc.opt.2[,2],
      lwd=2,col=col[3])
lines(ncc.mu_c.2[,1],ncc.mu_c.2[,2],
      lwd=2,col=col[4])
legend("bottomleft",
       legend=c("w/ CC (Climate Change)",
                "w/o CC + w/o optim.",
                "w/o CC + w/ optim.",
                "w/o CC + w/ optim. + imposed growth path"),
       lty=rep(1,4),
       col=col,
       lwd=rep(2,4),
       bty = "n",cex=2)
dev.off()

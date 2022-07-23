par(plt=c(.1,.9,.1,.9))

values.of.mu_d <- seq(0.4*model_sol$parameters$mu_d,
                      1.4*model_sol$parameters$mu_d,length.out = 8)
values.of.mu_n <- seq(0.4*model_sol$parameters$mu_n,
                      1.4*model_sol$parameters$mu_n,length.out = 8)


values.of.mu_d.mu_n <- cbind(values.of.mu_d %x% rep(1,length(values.of.mu_n)),
                             rep(1,length(values.of.mu_d)) %x% values.of.mu_n)


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
  
  all.Ts <- foreach(i = 1:dim(values.of.mu_d.mu_n)[1], .combine=rbind) %dopar% {
    
    model_i <- model_sol
    
    model_i$parameters$mu_d <- values.of.mu_d.mu_n[i,1]
    model_i$parameters$mu_n <- values.of.mu_d.mu_n[i,2]
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
    
    # SCC
    SCC.Q <- scc.fct(model_i_sol,0)
    SCC.P <- scc.fct(model_i_sol_P,0)
    
    c(T.P,T.Q,SCC.P,SCC.Q,T.P.check)
  }
  
  stopCluster(cl)
  
  
  T.P       <- matrix(all.Ts[,1],length(values.of.mu_n),length(values.of.mu_d))
  T.P.check <- matrix(all.Ts[,5],length(values.of.mu_n),length(values.of.mu_d))
  T.Q       <- matrix(all.Ts[,2],length(values.of.mu_n),length(values.of.mu_d))
  SCC.P     <- matrix(all.Ts[,3],length(values.of.mu_n),length(values.of.mu_d))
  SCC.Q     <- matrix(all.Ts[,4],length(values.of.mu_n),length(values.of.mu_d))
  
  RP.T   <- matrix((T.Q - T.P)/T.Q,length(values.of.mu_n),
                   length(values.of.mu_d))
  RP.SCC <- matrix((SCC.Q - SCC.P)/SCC.Q,length(values.of.mu_n),
                   length(values.of.mu_d))
  

#Plots
#1
FILE = paste("/outputs/Figures/Figure_RiskPremium.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=7, height=3)

par(mfrow=c(1,2))
par(plt=c(.15,.95,.15,.85))

# Refine grid (interpolation):
values.of.mu_d.fine <- seq(values.of.mu_d[1],
                           tail(values.of.mu_d,1),length.out = 40)
values.of.mu_n.fine <- seq(values.of.mu_n[1],
                           tail(values.of.mu_n,1),length.out = 40)

# Use splines:
x <- c(matrix(values.of.mu_n,length(values.of.mu_n),length(values.of.mu_d)))
y <- c(t(matrix(values.of.mu_d,length(values.of.mu_d),length(values.of.mu_n))))
mod <- mgcv::gam(c(T.P) ~ te(x, y))
T.P.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                  length(values.of.mu_n.fine),length(values.of.mu_d.fine))
mod <- mgcv::gam(c(T.Q) ~ te(x, y))
T.Q.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                  length(values.of.mu_n.fine),length(values.of.mu_d.fine))
mod <- mgcv::gam(c(RP.T) ~ te(x, y))
RP.T.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                   length(values.of.mu_n.fine),length(values.of.mu_d.fine))

mod <- mgcv::gam(c(SCC.P) ~ te(x, y))
SCC.P.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                    length(values.of.mu_n.fine),length(values.of.mu_d.fine))
mod <- mgcv::gam(c(SCC.Q) ~ te(x, y))
SCC.Q.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                    length(values.of.mu_n.fine),length(values.of.mu_d.fine))
mod <- mgcv::gam(c(RP.SCC) ~ te(x, y))
RP.SCC.fit <- matrix(predict(mod, expand.grid(x = values.of.mu_n.fine, y = values.of.mu_d.fine)),
                     length(values.of.mu_n.fine),length(values.of.mu_d.fine))


contour(values.of.mu_n.fine,values.of.mu_d.fine,T.P.fit,
        vfont = c("sans serif", "bold"),las=1,
        ylab=expression(paste("Disaster magnitude ",mu[D],sep="")),
        xlab=expression(paste("Feedback loops magnitude ",mu[N],sep="")),
        lwd=2,col=P.col.line,labcex=1,las=1,yaxt='n',
        main="(a) - Temperature expectations",
        cex.main=1.2,cex.axis=1,cex.lab=1.2
)
abline(h=model_sol$parameters$mu_d,col="grey",lty=3,lwd=2)
abline(v=model_sol$parameters$mu_n,col="grey",lty=3,lwd=2)
points(model_sol$parameters$mu_n, model_sol$parameters$mu_d,pch=15,
       col="dark grey",cex=1.6)
par(new=TRUE)
contour(values.of.mu_n.fine,values.of.mu_d.fine,T.Q.fit,
        vfont = c("sans serif", "bold"),
        ylab="",xlab="",las=1,
        lwd=2,col=Q.col.line,labcex=1)

contour(values.of.mu_n.fine,values.of.mu_d.fine,RP.T.fit,
        vfont = c("sans serif", "bold"),las=1,
        ylab=expression(paste("Disaster magnitude ",mu[D],sep="")),
        xlab=expression(paste("Feedback loops magnitude ",mu[N],sep="")),
        lwd=2,col="black",labcex=1,las=1,
        main="(b) - Share of risk premium in swap rate",
        cex.main=1.2,cex.axis=1,cex.lab=1.2)
abline(h=model_sol$parameters$mu_d,col="grey",lty=3,lwd=2)
abline(v=model_sol$parameters$mu_n,col="grey",lty=3,lwd=2)
points(model_sol$parameters$mu_n, model_sol$parameters$mu_d,pch=15,
       col="dark grey",cex=1.6)
dev.off()


#2
FILE = paste("/outputs/Figures/Figure_SCC.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=7, height=3)

par(mfrow=c(1,2))
par(plt=c(.15,.95,.15,.85))

contour(values.of.mu_n.fine,values.of.mu_d.fine,SCC.P.fit,
        vfont = c("sans serif", "bold"),
        ylab=expression(paste("Disaster magnitude ",mu[D],sep="")),
        xlab=expression(paste("Feedback loops magnitude ",mu[N],sep="")),
        lwd=2,col=P.col.line,labcex=1,lty=2,las=1,yaxt='n',
        main="(a) - Social Cost of Carbon (SCC)",
        cex.main=1.2,cex.axis=1,cex.lab=1.2)
abline(h=model_sol$parameters$mu_d,col="grey",lty=3,lwd=2)
abline(v=model_sol$parameters$mu_n,col="grey",lty=3,lwd=2)
points(model_sol$parameters$mu_n, model_sol$parameters$mu_d,pch=15,
       col="dark grey",cex=1.6)
par(new=TRUE)
contour(values.of.mu_n.fine,values.of.mu_d.fine,SCC.Q.fit,
        vfont = c("sans serif", "bold"),
        xlab="",ylab="",las=1,
        lwd=2,col=Q.col.line,labcex=1)

contour(values.of.mu_n.fine,values.of.mu_d.fine,RP.SCC.fit,
        vfont = c("sans serif", "bold"),
        ylab=expression(paste("Disaster magnitude ",mu[D],sep="")),
        xlab=expression(paste("Feedback loops magnitude ",mu[N],sep="")),
        lwd=2,col="black",labcex=1,las=1,
        main="(b) - Share of risk premium in SCC",
        cex.main=1.2,cex.axis=1,cex.lab=1.2)
abline(h=model_sol$parameters$mu_d,col="grey",lty=3,lwd=2)
abline(v=model_sol$parameters$mu_n,col="grey",lty=3,lwd=2)
points(model_sol$parameters$mu_n, model_sol$parameters$mu_d,pch=15,
       col="dark grey",cex=1.6)

dev.off()

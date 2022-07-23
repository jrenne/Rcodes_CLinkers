plot(ecdf(test$T_at[H,]),col="darkorange")

#RCP data
temp<-read.table("./data/mean_ssp.txt",header=TRUE)
temp_graph<-temp[3:11,]

#distribution and expected temperature
values.of.temperatures <- seq(as.integer(min(test$T_at[H,])),
                              ceiling(max(test$T_at[H,]))+1,
                              by=.2)

  cl <- makeCluster(number.of.cores)
  registerDoParallel(cl)
  
  save.image("outputs/parallel.Rdata")
  clusterEvalQ(cl,load("outputs/parallel.Rdata"))
  
  clusterEvalQ(cl,library(MASS))
  clusterEvalQ(cl,library(expm))
  
  print("Compute Q probas")
  all.Probas.Q <- foreach(i = 1:length(values.of.temperatures), 
                          .combine=rbind) %dopar% {
                            b    <- values.of.temperatures[i]
                            Temperature.Q<-c(varphi.hat(model_sol,
                                                        model_sol$omega_ZCB,
                                                        H,x,a,b)[[1]]/
                                               varphi(model_sol,model_sol$omega_ZCB,
                                                      H)[[3]])
                            Temperature.Q
                          }
  print("Compute P probas")
  
  
  all.Probas.P <- foreach(i = 1:H, .combine=rbind) %dopar% {
    Temperature.P <- fourier(model_sol,x,values.of.temperatures,i,9)
    Temperature.P
  }
  
  stopCluster(cl)
  
  all.Probas.P<-t(all.Probas.P)
  
  nb.temperatures <- 200                                                        #number of temp levels
  scale.temperatures.values <- seq(values.of.temperatures[1],
                                   tail(values.of.temperatures,1),
                                   length.out = nb.temperatures+1)
  all.pdf.Q <- NULL
  all.cdf.Q <- NULL
  all.pdf.P <- NULL
  all.cdf.P <- NULL
  
  for(t in 1:H){
    # make sure probas are increasing:
    cdf.Q <- pmax(pmin(all.Probas.Q[,t],1),0)
    cdf.P <- pmax(pmin(all.Probas.P[,t],1),0)
    for(i in 2:length(cdf.P)){
      if(cdf.Q[i]<cdf.Q[i-1]){
        cdf.Q[i] <- cdf.Q[i-1]
      }
      if(cdf.P[i]<cdf.P[i-1]){
        cdf.P[i] <- cdf.P[i-1]
      }
    }
    
    tmp.Q <- splinefun(x=values.of.temperatures, y=cdf.Q, method="hyman")
    tmp.P <- splinefun(x=values.of.temperatures, y=cdf.P, method="hyman")
    
    fitted.cdf.values.Q <- tmp.Q(scale.temperatures.values)
    fitted.pdf.values.Q <- diff(fitted.cdf.values.Q)
    fitted.cdf.values.P <- tmp.P(scale.temperatures.values)
    fitted.pdf.values.P <- diff(fitted.cdf.values.P)
    
    all.pdf.Q <- cbind(all.pdf.Q,fitted.pdf.values.Q)
    all.cdf.Q <- cbind(all.cdf.Q,fitted.cdf.values.Q)
    
    all.pdf.P <- cbind(all.pdf.P,fitted.pdf.values.P)
    all.cdf.P <- cbind(all.cdf.P,fitted.cdf.values.P)
  }
  
  # Computation of confidence intervals:
  get.index.CI <- function(q,cdf){
    lower.q <- (1-q)/2
    upper.q <- q + (1-q)/2
    # get index for lower q:
    lower.index <- which(cdf>lower.q)[1]
    upper.index <- which(cdf>upper.q)[1]
    return(c(lower.index,upper.index))
  }
  all.CI.Q <- array(NaN,c(2,dim(all.cdf.Q)[2],length(vector.of.CI)))
  all.CI.P <- array(NaN,c(2,dim(all.cdf.P)[2],length(vector.of.CI)))
  for(i in 1:length(vector.of.CI)){
    q <- vector.of.CI[i]
    for(j in 1:dim(all.cdf.P)[2]){
      aux <- get.index.CI(q,all.cdf.Q[,j])
      all.CI.Q[1,j,i] <- aux[1] # lower bound
      all.CI.Q[2,j,i] <- aux[2] # upper bound
      aux <- get.index.CI(q,all.cdf.P[,j])
      all.CI.P[1,j,i] <- aux[1] # lower bound
      all.CI.P[2,j,i] <- aux[2] # upper bound
    }
  }
  Tplot<-list(all.cdf.P,all.cdf.Q,all.CI.P,all.CI.Q,all.pdf.P,all.pdf.Q,
              scale.temperatures.values,nb.temperatures)
  save(Tplot,file="outputs/backup_parallel/pdf_temp.Rdat")

# Compute expected trajectories (under P and Q):
EV      <- EV.fct(model_sol,h=H)
ET.P    <- EV$EX$T_at[1:H]
ET.Q    <- varphi.tilde(model_sol,omega_T.at,H)[[1]]/
  varphi(model_sol,omega_ZCB,H)[[3]]

#Plots
FILE = paste("/outputs/Figures/Figure_Tat_P_and_Q_vector_CI.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=9, height=6)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1,2), heights=c(1,1,2))
par(plt=c(.15,.95,.15,.85))

y.lim <- c(1.5,6)
x.lim <- c(2040,2100)

plot(model_sol$vec_date[2:(H+1)],ET.P,
     xlim=x.lim,ylim=y.lim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",xlab="",ylab="Degrees Celsius",las=1,
     main="(a) - Trajectory of atm. temperature")

for(i in length(vector.of.CI):1){
  # P
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          scale.temperatures.values[c(all.CI.P[1,,i],rev(all.CI.P[2,,i]))],
          col=P.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      ET.P,lwd=2,col=P.col.line)
lines(temp_graph[,1],temp_graph[,5],lty=2,lwd=1,col="lightsteelblue4")
lines(temp_graph[,1],temp_graph[,4],lty=2,lwd=1,col="lightsteelblue4")
lines(temp_graph[,1],temp_graph[,2],lty=2,lwd=2,col="grey28")
legend("topright",
       legend=c("RCP4.5+RCP6.0","+/- 2 std"),
       lty=c(2,2),
       col=c("grey28","lightsteelblue4"),cex=1.5,
       lwd=c(2,1),bty = "n")

plot(model_sol$vec_date[2:(H+1)],ET.Q,
     xlim=x.lim,ylim=y.lim,las=1,
     col="white",xlab="",ylab="Degrees Celsius", cex.main=1.5,cex.axis=1.5,
     cex.lab=1.5,
     main="(b) - Trajectory of atm. temperature including Risk-Premium")
for(i in length(vector.of.CI):1){
  # Q
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          scale.temperatures.values[c(all.CI.Q[1,,i],rev(all.CI.Q[2,,i]))],
          col=Q.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      ET.Q,lwd=2,col=Q.col.line)
lines(model_sol$vec_date[2:(H+1)],
      ET.P,lwd=2,col=P.col.line)

plot(scale.temperatures.values[2:(nb.temperatures+1)],all.pdf.P[,H],type="l",
     col=P.col.line,lwd=3,xlim=c(2,max(scale.temperatures.values)),
     main="(c) - P.d.f. of atm. temperature in 2100",cex.main=1.5,cex.axis=1.5,
     cex.lab=1.5,
     xlab="",ylab="",yaxt="no")
lines(scale.temperatures.values[2:(nb.temperatures+1)],all.pdf.Q[,H],
      col=Q.col.line,lwd=3)
abline(v=scale.temperatures.values[which.min(abs(all.cdf.Q[,model_sol$horiz.2100]
                                                 -0.5))],
       col=Q.col.line,lty=3,lwd=2)                                              #median of Q
abline(v=ET.Q[H],col=Q.col.line,lty=1,lwd=1)                                    #mean of Q
abline(v=scale.temperatures.values[which.min(abs(all.cdf.P[,model_sol$horiz.2100]
                                                 -0.5))],
       col=P.col.line,lty=3,lwd=2)                                              #median of P
abline(v=ET.P[H],col=P.col.line,lty=1,lwd=1)                                    #mean of P

legend("topright",
       legend=c("Physical p.d.f.","Risk-Adjusted p.d.f.","Mean","Median"),
       lty=c(1,1,1,3),
       col=c(P.col.line,Q.col.line,"black","black"),cex=1.5,
       lwd=c(2,2,1,2),bty = "n")
dev.off()

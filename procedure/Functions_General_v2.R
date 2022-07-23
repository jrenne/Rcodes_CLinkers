#Generalized Functions file
#ctrl F: CHANGE
################################################################################*MISCELLANEOUS
###Plot with confidence area

make_confidence_area <- function(x,y,pdf_xy,p,tol=10^(-2)){
  # Computes the coordinates of a polygon that delineates a confidence
  #     area associated with probability p.
  # The components of x (and of y) must be equally spaced.
  
  h.x <- x[2] - x[1]
  h.y <- y[2] - y[1]
  
  pdf_xy_h2 <- pdf_xy * h.x * h.y # The sum of all entries of pdf_xy_h2 should ne close to 1.
  
  # Use bisection to find limit value of pdf_xy_h2 to get probability p within the area:
  z.up <- max(pdf_xy_h2)
  z.lo <- 0
  error <- 1000
  while(error > tol){
    pdf_xy_h2_aux <- pdf_xy_h2
    pdf_xy_h2_aux[pdf_xy_h2<(z.lo+z.up)/2] <- 0
    aux <- sum(pdf_xy_h2_aux) - p
    if(aux>0){
      z.lo <- (z.lo+z.up)/2
    }else{
      z.up <- (z.lo+z.up)/2
    }
    error <- abs(aux)
    print(error)
  }
  
  M <- pdf_xy_h2_aux
  M[M>0] <- 1 # Then the matrix M is filled only with 0 and 1.
  x.polygon.1 <- NULL
  x.polygon.2 <- NULL
  y.polygon.1 <- NULL
  y.polygon.2 <- NULL
  for(i in 1:length(x)){
    if(sum(M[i,])>0){
      x.polygon.1 <- c(x.polygon.1,x[i])
      x.polygon.2 <- c(x.polygon.2,x[i])
      index.y <- which(M[i,]>0)
      y.polygon.1 <- c(y.polygon.1,y[index.y[1]])
      y.polygon.2 <- c(y.polygon.2,y[index.y[length(index.y)]])
    }
  }
  
  x.polygon <- c(x.polygon.1,rev(x.polygon.2),x.polygon.1[1])
  y.polygon <- c(y.polygon.1,rev(y.polygon.2),y.polygon.1[1])
  
  return(list(
    x.polygon = x.polygon,
    y.polygon = y.polygon
  ))
}

#######################################
#Function to extract n'th element of the list

extract<-function(list, n){
  sapply(list, `[`, n)
}

################################################################################*STOCKS
##Activation Cum_div
#* gives initial value for mu_pd
div.fct<-function(model_sol,h=NaN,omega.PDh=model_sol_s$omega_Div,
                  theta=theta0,indic_mitig=FALSE,X=model_sol_s$X,t=0,
                  mu.chosen=model_sol_s$mu){
  EV                <-EV.fct(model_sol,h)
  param             <-model_sol$parameters
  model_sol$mu_div.0<-matrix(c((1-param$chi.div)*X[1],
                               (1-param$chi.div)*EV$EX$delc),
                             length(EV$EX$delc)+1,1)
  
  model_sol_div     <- model_solve(model_sol,theta,indic_mitig,mu.chosen)
  if(is.nan(h)){
    PD_h            <- varphi(model_sol_div,omega.PDh,h-2,X,t)
  }else{
    PD_h            <- varphi(model_sol_div,omega.PDh,h,X,t)$P.t
  }
  
  PD                <- sum(PD_h)
  pd                <- log(PD)
  
  model_sol_div[["PD_h"]]<-PD_h
  model_sol_div[["PD"]]  <-PD
  model_sol_div[["pd"]]  <-pd
  
  model_sol_div[["mu_div.1"]]<-matrix(c(param$chi.div,
                                        rep(0,model_sol$n.X-1)),
                                      model_sol$n.X,1)
  
  model_sol_div$ini_matx[["mu_div.0"]]<-model_sol_div$mu_div.0[1]
  model_sol_div$mu_div.0              <-model_sol_div$mu_div.0[-1]
  
  return(model_sol_div)
}

#returns a column vector with kappa0,kappa1
kappa.fct<-function(model_sol,pd_bar){
  kappa1<-exp(pd_bar)/(1+exp(pd_bar))
  kappa0<-log(1+exp(pd_bar))-kappa1*pd_bar
  kappa <-rbind(kappa0,kappa1)
  
  return(kappa)
}

#mu_pd function for GN algorithm/optim
#*mu_pd is a 22*1 matrix
mu_pd4GN.fct<-function(model_sol,mu_pd,EV){
  
  EX.max <-matrix(extract(EV$EX,length(EV$date)),ncol=1) #changer length avec Tmax/ne pas prévoir plus loin que Tmax
  VX.max <-EV$CovX[[length(EV$date)]]
  
  pd_bar  <-mu_pd[1]+t(mu_pd[2:(model_sol$n.X+1)])%*%EX.max+
    1/2*t(mu_pd[2:(model_sol$n.X+1)])%*%VX.max%*%mu_pd[2:(model_sol$n.X+1)] 
  
  kap     <-kappa.fct(model_sol,pd_bar)
  
  mu_pd1.right<--model_sol$inf_matx$eta1.inf+
    b1.fct.inf(model_sol,
               kap[2]*mu_pd[2:(model_sol$n.X+1)]+
                 model_sol$mu_div.1+
                 model_sol$inf_matx$pi.inf)-
    b1.fct.inf(model_sol,model_sol$inf_matx$pi.inf)
  
  mu_pd0.right<-1/(1-kap[2])*
    (kap[1]+model_sol$mu_div.0[length(model_sol$mu_div.0)]-
       model_sol$inf_matx$eta0.inf+
       a1.fct.inf(model_sol,
                  kap[2]*mu_pd[2:(length(model_sol$X)+1)]+
                    model_sol$mu_div.1+model_sol$inf_matx$pi.inf)-
       a1.fct.inf(model_sol,model_sol$inf_matx$pi.inf))
  
  
  mu_pd.right <-rbind(mu_pd0.right, mu_pd1.right)
  dev.mu_pd   <-mu_pd-mu_pd.right
  
  return(list("dev"=dev.mu_pd,"pd_bar"=pd_bar,"kappa"=kap))
}

######
lossfunction4stock <- function(model,mu_pd.complete,mu_pd.reduced,EV,FILTER){
  
  mu_pd           <-mu_pd.complete
  mu_pd[FILTER==1]<-mu_pd.reduced
  
  dev<-mu_pd4GN.fct(model,mu_pd,EV)$dev[FILTER==1,,drop=FALSE]
  
  
  # Compute loss function:
  loss <- sum(dev^2)*10^(4)
  return(loss)
}

######
stock.solve<-function(model_sol,h,omega.PDh=model_sol$omega_Div,
                      theta=theta0,X=model_sol$X,t=0,
                      mu.chosen=model_sol$mu,indic_mitig=FALSE){
  ##################
  #Data preparation#
  ##################
  #state vector nb var                                                            #CHANGE
  
  #model_sol$parameters$chi.div<-3 #!!!!
  model_sol_s     <-model_sol
  model_sol_s$n.Z <-model_sol$n.Z+5 #Cum_div,pd,r.s,Cum_rs,r.1
  
  model_sol_s$n.W <-model_sol$n.W+1 #Div
  
  #Update vector elements in param
  filter<-list("D"=c(which(model_sol_s$parameters$ell1.D!=0),
                     model_sol_s$parameters$ell1.D[model_sol_s$parameters$ell1.D!=0]),
               "N"=c(which(model_sol_s$parameters$ell1.N!=0),
                     model_sol_s$parameters$ell1.N[model_sol_s$parameters$ell1.N!=0]),
               "div"=c(which(model_sol_s$parameters$ell1.div!=0),
                       model_sol_s$parameters$ell1.div[model_sol_s$parameters$ell1.div!=0]))
  
  model_sol_s$parameters$ell1.D  <-matrix(0,nrow=model_sol_s$n.Z+model_sol_s$n.W)
  model_sol_s$parameters$ell1.N  <-matrix(0,nrow=model_sol_s$n.Z+model_sol_s$n.W)
  model_sol_s$parameters$ell1.div<-matrix(0,nrow=model_sol_s$n.Z+model_sol_s$n.W)
  
  model_sol_s$parameters$ell1.D[filter$D[1]]    <-filter$D[2]
  model_sol_s$parameters$ell1.N[filter$N[1]]    <-filter$N[2]
  model_sol_s$parameters$ell1.div[filter$div[1]]<-filter$div[2]
  
  model_sol_s$ell1.N.tilde <-model_sol_s$parameters$ell1.N
  model_sol_s$ell1.N.tilde[model_sol_s$n.Z+model_sol_s$n.eta+2]<-param$rho.N/param$mu_n
  
  #Update vector elements in load_ini_model
  model_sol_s$mu_r1.a1<-rep(list(matrix(0,nrow=model_sol_s$n.Z+model_sol_s$n.W)),
                            model_sol_s$Tmax)
  model_sol_s$mu_pd.a0<-rep(list(matrix(0,nrow=model_sol_s$n.Z+model_sol_s$n.W+1)),
                            model_sol_s$Tmax)
  
  #Update pricing vectors
  model_sol_s$omega_ZCB     <-matrix(0,model_sol_s$n.Z+model_sol_s$n.W,1)
  model_sol_s$omega_T.at    <-model_sol_s$omega_ZCB
  model_sol_s$omega_T.at[9] <-1
  model_sol_s$omega_Div     <-model_sol_s$omega_ZCB
  model_sol_s$omega_Div[15] <-1
  
  #Update model_sol w/ extended matrices
  model_sol_tempo<-model_solve(model_sol_s,theta,
                               indic_mitig,mu.chosen)
  
  model_sol<-model_sol_tempo
  #Solve the model w/ dividends (to construct pd)
  model_sol_stock<-model_sol
  model_sol_stock<-div.fct(model_sol_stock,max(h,model_sol_stock$Tmax),
                           omega.PDh,theta,indic_mitig=FALSE,X,t,mu.chosen)
  ###
  ##construction of mu_pd for pd
  ###
  #Time indep extension
  #####
  if(h>=(model_sol$Tmax-1)){
    pi      <-model_sol_stock$pi
    inf     <-rep(list(model_sol_stock$inf_matx$pi.inf)
                  ,h-(model_sol_stock$Tmax-1)+1)
    pi      <-c(pi,inf)
    
    eta_1   <-model_sol_stock$eta1
    inf     <-rep(list(model_sol_stock$inf_matx$eta1.inf),
                  h-(model_sol_stock$Tmax-1)+1)
    eta_1   <-c(eta_1,inf)
    
    eta_0   <-model_sol_stock$eta0
    inf     <-matrix(model_sol_stock$inf_matx$eta0.inf,
                     h-(model_sol_stock$Tmax-1)+1,1)
    eta_0   <-rbind(eta_0,inf)
  }else{
    pi      <-model_sol_stock$pi
    eta_1   <-model_sol_stock$eta1
    eta_0   <-model_sol_stock$eta0
  }
  #####
  #optimize only the relevant elements of the state vector                      #CHANGE if add persistent shock
  FILTER2        <-c(1,rep(0,model_sol_stock$n.X))                              #1 is for mu_pd0
  FILTER2[4]     <-1 #E
  FILTER2[7:11]  <-1 #M,T
  
  FILTER2[model_sol_stock$n.Z+1]                       <-1 #etaF
  FILTER2[model_sol_stock$n.Z+model_sol_stock$n.eta+2] <-1 #N
  
  #Update expectations and variance after adding Cum_div
  EV<-EV.fct(model_sol_stock,max(h,model_sol_stock$Tmax))
  
  #create an approximated initial values, by linearizing our formulas
  list_ini.GN            <-list("pd"=model_sol_stock$pd)
  list_ini.GN[["kap"]]   <-kappa.fct(model_sol_stock,list_ini.GN$pd)
  list_ini.GN[["mu_pd1"]]<-solve(diag(model_sol_stock$n.X)-
                                   list_ini.GN$kap[2]*model_sol_stock$M)%*%
    (-model_sol_stock$inf_matx$eta1.inf+model_sol_stock$M%*%model_sol_stock$mu_div.1)
  
  
  model_sol_stock[["list_ini.GN.div"]]<-list_ini.GN
  
  #construct mu_pd=[mu_pd0,mu_pd1] with our initial guess
  mu_pd.right   <-rbind(1/(1-list_ini.GN$kap[2])*
                          (list_ini.GN$kap[1]+
                             model_sol_stock$mu_div.0[length(model_sol_stock$mu_div.0)]-
                             model_sol_stock$inf_matx$eta0.inf+
                             a1.fct.inf(model_sol_stock,
                                        list_ini.GN$kap[2]*list_ini.GN$mu_pd1+
                                          model_sol_stock$mu_div.1+
                                          model_sol_stock$inf_matx$pi.inf)-
                             a1.fct.inf(model_sol_stock,
                                        model_sol_stock$inf_matx$pi.inf)),
                        list_ini.GN$mu_pd1)
  #prepare for optimization
  mu_pd.inf  <-mu_pd.right
  for(i in 1:7){
    STOCK.optim<-optimx(c(mu_pd.inf[FILTER2==1]),
                        lossfunction4stock,
                        model = model_sol_stock,
                        EV=EV,
                        mu_pd.complete=mu_pd.inf,
                        FILTER=FILTER2,
                        method = "nlminb",
                        control=list(trace=1,
                                     kkt = FALSE))
    mu_pd.filter         <-c(as.matrix(STOCK.optim)[1:sum(FILTER2)])
    mu_pd.inf            <-matrix(0,model_sol_stock$n.X+1)
    mu_pd.inf[FILTER2==1]<-mu_pd.filter
    STOCK.optim<-optimx(c(mu_pd.inf[FILTER2==1]),
                        lossfunction4stock,
                        model = model_sol_stock,
                        EV=EV,
                        mu_pd.complete=mu_pd.inf,
                        FILTER=FILTER2,
                        method = "Nelder-Mead",
                        control=list(trace=1,maxit=10000,
                                     kkt = FALSE))
    mu_pd.filter         <-c(as.matrix(STOCK.optim)[1:sum(FILTER2)])
    mu_pd.inf            <-matrix(0,model_sol_stock$n.X+1)
    mu_pd.inf[FILTER2==1]<-mu_pd.filter
  }
  #better initial guess, restart optimiz w/o the approx
  #####
  # STOCK.optim<-optimx(c(mu_pd.right[FILTER2==1]),
  #                     lossfunction4stock,
  #                     model = model_sol_stock,
  #                     EV=EV,
  #                     mu_pd.complete=mu_pd.right,
  #                     FILTER=FILTER2,
  #                     method = "nlminb",
  #                     control=list(trace=1,
  #                                  kkt = FALSE))
  # mu_pd.filter         <-c(as.matrix(STOCK.optim)[1:sum(FILTER2)])
  # mu_pd.inf            <-matrix(0,model_sol_stock$n.X+1)
  # mu_pd.inf[FILTER2==1]<-mu_pd.filter
  #####
  #recursive mu_pd
  mu_pd<-list()
  pd.b <-list()
  kappa<-list()
  #start values
  model_sol_stock$ini_matx[["pd.b"]] <-model_sol_stock$pd
  model_sol_stock$ini_matx[["kappa"]]<-kappa.fct(model_sol_stock,model_sol_stock$pd)
  #last values
  mu_pd[[h]]  <-mu_pd.inf
  pd.b[[h]]   <-mu_pd4GN.fct(model_sol_stock,mu_pd[[h]],EV)$pd_bar
  kappa[[h]]  <-mu_pd4GN.fct(model_sol_stock,mu_pd[[h]],EV)$kappa
  kap         <-kappa[[h]]
  
  #case where h-i>98, need time-indep formulas
  i<-1
  while((h-i)>=(model_sol_stock$Tmax-1)){
    EXh     <-matrix(extract(EV$EX,h-i),ncol=1)
    VXh     <-EV$CovX[[h-i]]
    
    mu_pd[[h-i]]<-rbind(kap[1]+kap[2]*mu_pd[[h-i+1]][1]+
                          model_sol_stock$mu_div.0[h-i+1]-eta_0[h-i+1]+
                          a1.fct.inf(model_sol_stock,pi[[h-i+1]]+kap[2]*
                                       mu_pd[[h-i+1]][2:(model_sol_stock$n.X+1)]+
                                       model_sol_stock$mu_div.1)-
                          a1.fct.inf(model_sol_stock,pi[[h-i+1]]),
                        -eta_1[[h-i+1]]+
                          b1.fct.inf(model_sol_stock,pi[[h-i+1]]+kap[2]*
                                       mu_pd[[h-i+1]][2:(model_sol_stock$n.X+1)]+
                                       model_sol_stock$mu_div.1)-
                          b1.fct.inf(model_sol_stock,pi[[h-i+1]]))
    
    pd_bar  <-mu_pd[[h-i]][1]+
      t(mu_pd[[h-i]][2:(model_sol_stock$n.X+1)])%*%EXh+
      1/2*t(mu_pd[[h-i]][2:(model_sol_stock$n.X+1)])%*%
      VXh%*%mu_pd[[h-i]][2:(model_sol_stock$n.X+1)]
    
    pd.b[[h-i]]<-pd_bar
    
    kap         <-kappa.fct(model_sol_stock,pd_bar)
    kappa[[h-i]]<-c(kap)
    
    i<-i+1
  }
  #h-k<98, time-dep formulas
  for(k in i:(h-1)){
    mu_pd[[h-k]]<-rbind(kap[1]+kap[2]*mu_pd[[h-k+1]][1]+
                          model_sol_stock$mu_div.0[h-k+1]-
                          eta_0[h-k+1]+
                          a1.fct.1(model_sol_stock,pi[[h-k+1]]+kap[2]*
                                     mu_pd[[h-k+1]][2:(model_sol_stock$n.X+1)]+
                                     model_sol_stock$mu_div.1,h-k)-
                          a1.fct.1(model_sol_stock,pi[[h-k+1]],h-k),
                        -eta_1[[h-k+1]]+
                          b1.fct.1(model_sol_stock,pi[[h-k+1]]+kap[2]*
                                     mu_pd[[h-k+1]][2:(model_sol_stock$n.X+1)]+
                                     model_sol_stock$mu_div.1,
                                   h-k)-
                          b1.fct.1(model_sol_stock,pi[[h-k]],h-k))
    
    EXh     <-matrix(extract(EV$EX,h-k),ncol=1)
    VXh     <-EV$CovX[[h-k]]
    
    pd_bar     <-mu_pd[[h-k]][1]+
      t(mu_pd[[h-k]][2:(model_sol_stock$n.X+1)])%*%EXh+
      1/2*t(mu_pd[[h-k]][2:(model_sol_stock$n.X+1)])%*%
      VXh%*%mu_pd[[h-k]][2:(model_sol_stock$n.X+1)]
    
    pd.b[[h-k]]<-pd_bar
    
    kap         <-kappa.fct(model_sol_stock,pd_bar)
    kappa[[h-k]]<-c(kap)
  }
  
  model_sol_stock[["mu_pd"]]         <-mu_pd
  model_sol_stock[["pd.b.full"]]     <-matrix(unlist(pd.b),ncol=1)
  model_sol_stock[["kap.full"]]      <-kappa
  model_sol_stock$ini_matx[["mu_pd"]]<-rbind(kappa[[1]][1]+kappa[[1]][2]*mu_pd[[1]][1]+
                                               model_sol_stock$mu_div.0[1]-
                                               eta_0[1]+
                                               a1.fct.1(model_sol_stock,pi[[1]]+kappa[[1]][2]*
                                                          mu_pd[[1]][2:(model_sol_stock$n.X+1)]+
                                                          model_sol_stock$mu_div.1,0)-
                                               a1.fct.1(model_sol_stock,pi[[1]],0),
                                             -eta_1[[1]]+
                                               b1.fct.1(model_sol_stock,pi[[1]]+kappa[[1]][2]*
                                                          mu_pd[[1]][2:(model_sol_stock$n.X+1)]+
                                                          model_sol_stock$mu_div.1,
                                                        0)-
                                               b1.fct.1(model_sol_stock,pi[[1]],0))
  
  model_sol_stock$ini_matx[["pd"]]   <-model_sol_stock$ini_matx$mu_pd[1]+
    t(model_sol_stock$ini_matx$mu_pd[2:(model_sol_stock$n.X+1)])%*%X
  model_sol_stock$ini_pd             <-model_sol_stock$ini_matx$pd
  model_sol_stock$ini_matx[["kap"]]  <-kappa.fct(model_sol_stock,
                                                 model_sol_stock$ini_matx$pd)
  
  #construct stock returns
  model_sol_stock$mu_div.0 <-c(model_sol_stock$ini_matx$mu_div.0,
                               model_sol_stock$mu_div.0)
  model_sol_stock$kap0.om0 <-c(model_sol_stock$ini_matx$kap[1],
                               extract(model_sol_stock$kap.full,1))
  model_sol_stock$kap1.a0  <-c(model_sol_stock$ini_matx$kap[2],
                               extract(model_sol_stock$kap.full,2))
  model_sol_stock$mu_pd.a0 <-c(list(model_sol_stock$ini_matx$mu_pd),
                               model_sol_stock$mu_pd)
  model_sol_stock$mu_r0.om0<-c(0,
                               model_sol_stock$eta0,
                               model_sol_stock$inf_matx$eta0.inf)
  model_sol_stock$mu_r1.a1 <-c(list(matrix(0,model_sol_stock$n.X,1)),
                               model_sol_stock$eta1,
                               list(model_sol_stock$inf_matx$eta1.inf))
  if(h<model_sol_stock$Tmax){
    model_sol_stock$kap0.om0 <-c(model_sol_stock$kap0.om0,
                                 rep(model_sol_stock$kap0.om0[h],
                                     model_sol_stock$Tmax-h))
    model_sol_stock$kap1.a0  <-c(model_sol_stock$kap1.a0,
                                 rep(model_sol_stock$kap1.a0[h],
                                     model_sol_stock$Tmax-h))
    model_sol_stock$mu_pd.a0 <-c(model_sol_stock$mu_pd.a0,
                                 rep(model_sol_stock$mu_pd.a0[[h]],
                                     model_sol_stock$Tmax-h))
    model_sol_stock$mu_r0.om0<-c(model_sol_stock$mu_r0.om0,
                                 rep(model_sol_stock$mu_r0.om0[h],
                                     model_sol_stock$Tmax-h))
    model_sol_stock$mu_r1.a1 <-c(model_sol_stock$mu_r1.a1,
                                 rep(model_sol_stock$mu_r1.a1[[h]],
                                     model_sol_stock$Tmax-h))
  }
  
  model_stock<-model_solve(model_sol_stock,theta,indic_mitig=FALSE,mu.chosen)
  
  #Delete first element t=0
  model_stock$mu_div.0 <-model_stock$mu_div.0[-1]
  model_stock$kap0.om0 <-model_stock$kap0.om0[-1]
  model_stock$kap1.a0  <-model_stock$kap1.a0[-1]
  model_stock$mu_pd.a0 <-model_stock$mu_pd.a0[-1]
  model_stock$mu_r0.om0<-model_stock$mu_r0.om0[2:(Tmax-1)]
  model_stock$mu_r1.a1 <-model_stock$mu_r1.a1[2:(Tmax-1)]
  
  return(model_stock)
  
}
##Simulations (nb.traj) of real and approximated pd. 
#*Returns graphical illustrations, as well as the avg of the simulations.
pd.real.fct<-function(model_sol,h,nb.traj,
                      omega.PDh=model_sol$omega_Div,X=model_sol$X,t=0){
  
  X.simul<-simul.function(model_sol,h+2*h,nb.traj)$X
  tic("Real pd computations")
  PD_h<- c(list(matrix(
    rep(varphi(model_sol,omega.PDh,h,X,0)$P.t,nb.traj),h,nb.traj)),
    lapply(1:(h-1),function(t){
      apply(X.simul[[t]],2,function(x){
        varphi(model_sol,
               omega.PDh,
               h,
               x,t)$P.t/exp(x[15])
      })
    }
    ))
  
  PD                <- matrix(unlist(lapply(1:length(PD_h),function(x){
    apply(PD_h[[x]],2,sum)
  })),h,nb.traj,byrow=TRUE)
  pd.real           <- log(PD)
  toc()
  pd.approx         <-matrix(unlist(lapply(1:(h-1),function(t){
    apply(X.simul[[t]],2,function(x){
      model_sol$mu_pd[[t]][1]+
        t(model_sol$mu_pd[[t]][2:(model_sol$n.X+1)])%*%x
    })
  }
  ))
  ,ncol=nb.traj,byrow=T)
  
  model_sol[["pd.T"]]<-pd.real
  model_sol[["pd.A"]]<-pd.approx
  model_sol$ini_matx[["pd.A"]]<-model_sol$ini_matx$mu_pd[1]+
    t(model_sol$ini_matx$mu_pd[2:(model_sol$n.X+1)])%*%X
  
  col<-plasma(nb.traj)
  #1
  par(mfrow=c(2,ceiling(nb.traj/2)))
  for(i in 1:nb.traj){
    plot(model_sol$pd.b.full,type="l",
         ylim=range(c(model_sol$pd.b.full,model_sol$pd.T,model_sol$pd.A)))
    lines(c(model_sol$ini_matx$pd.A,model_sol$pd.A[,i]),col="red")
    lines(model_sol$pd.T[,i],col=col[i]) 
  }
  #2
  par(mfrow=c(1,1))
  plot(apply(model_sol$pd.T,1,mean)-
         c(model_sol$ini_matx$pd.A,apply(model_sol$pd.A,1,mean)),
       type="l",col="red")
  abline(h=0,col="black",lwd=2)
  #3
  plot(apply(model_sol$pd.T,1,mean)-
         c(model_sol$ini_matx$pd.A,apply(model_sol$pd.A,1,mean)),
       type="l",col="white")
  abline(h=0,col="black",lwd=2)
  for(i in 1:nb.traj){
    lines(model_sol$pd.T[,i]-c(model_sol$ini_matx$pd.A,model_sol$pd.A[,i]),
          col=col[i])
  }
  #4
  plot(model_sol$pd.b.full,type="l",
       ylim=range(c(model_sol$pd.b.full,model_sol$pd.T,model_sol$pd.A)))
  lines(c(model_sol$ini_matx$pd.A,model_sol$pd.A),col="red")
  lines(apply(model_sol$pd.T,1,mean),col="darkorange1")
  
  return(model_sol)
}

################################################################################*SCC
##SCC, for h>=0
#*if all = TRUE, vector for h=0 to h.

scc.fct<-function(model_sol,h,all=F,mat=6,C_0=model_sol$c0,X=model_sol$X,t=0){    
  mu.u1.c    <-abs(extract(model_sol$mu_u1.t1,mat))
  if(h>(model_sol$Tmax-1)){
    mu.u1.c    <-c(mu.u1.c,rep(abs(model_sol$inf_matx$mu_u1[mat]),
                               h-model_sol$Tmax+1)) 
  }
  omega_c    <-matrix(0,model_sol$n.X,1)
  omega_c[13]<-1
  
  if(h==0){
    scc<--C_0*mu_u.t.fct(model_sol)[[2]][mat]*10^(3)
  } else{
    scc<-multi.lt.fct(model_sol,omega_c,h,X,t)$uX_t.h*
                                            mu.u1.c[h+1]*C_0*10^(3)
  }
  
  if(all){
    scc.all<--C_0*mu_u.t.fct(model_sol)[[2]][mat]*10^(3)
    for(i in 1:(h-1)){
      scc.all[i+1]<-multi.lt.fct(model_sol,omega_c,i,X)$uX_t.h*
        mu.u1.c[i+1]*C_0*10^(3)
    }
    return(scc.all)/(1-model_sol$parameters$delta)
  }
  return(scc/(1-model_sol$parameters$delta))
}

################################################################################*PRICING
## Stock returns with new model and expected returns over time h

compute_rs<-function(model_sol,h,X,t=0,mu_div.0=NaN,mu_div.1=NaN){
  if(!is.na(mu_div.1[1])){
    model_sol$mu_div.0<-mu_div.0
    model_sol$mu_div.1<-mu_div.1
  }
  model_stock<-stock.solve(model_sol,h,X=X,t=t)
  EV.rs      <-EV.fct(model_sol,h)$EX$r.s
  
  mylist<-list(model_stock,EV.rs)
  
  return(mylist)
}

## Yield curve for ZCB
#* Compute yield curve for a ZCB of maturity h
#* returns a (h*2) matrix, 
#* where the first column corresponds to the date,
#* and the second column corresponds to the rate in percent.

compute_yc<-function(model_sol,h,X=model_sol$X,t=0){
  yds  <-varphi(model_sol,omega.varphi = model_sol$omega_ZCB,h,X,t)
  yds.r<-yds$r.t
  yds.d<-yds$date
  
  r<-cbind(yds.d,yds.r)
  
  return(r)
}

## Cst maturity for ZCB
#* Compute constant maturity 'h' for ZCB 'nb' times.
#* returns a (h*2) matrix, 
#* where the first column corresponds to the date,
#* and the second column corresponds to the rate in percent.

compute_cst_h<-function(model_sol,h,nb,X=model_sol$X,t=0){
  EX<-cbind(X,t(extract(EV.fct(model_sol,nb)$EX,(t+1):nb)))
  r<-0
  for(j in t:(nb-1)){
    yield <-varphi(model_sol,omega.varphi = model_sol$omega_ZCB,h,EX[,j+1],j)
    r[j+1]<-yield$r.t[h]
  }
  date<-seq(model_sol$tstep*t+model_sol$vec_date[1],
            by=model_sol$tstep,
            length=nb)
  
  dr<-cbind(date,r)
  
  return(dr)
}

##Corollary 4: TIBs-leverage function
#*i is the variable indexed to the bond
#*H is the maturity
#*T0 is the indexed decided by the issuer at date t
#*chi is the leverage factor

TIB<-function(model_sol,chi,T0,H,i,X=model_sol$X,t=0){
  o.ZCB<-matrix(0,model_sol$n.X,1)
  o.t   <-o.ZCB
  o.t[i]<-1
  P.tib<-(1-chi*T0)*varphi(model_sol,o.ZCB,H,X,t)[[3]]+
          chi*varphi.tilde(model_sol,o.t,H,X,t)[[1]]
  r.tib<--log(P.tib)/(model_sol$tstep*((t+1):(t+H)))*100
  
  vec_date<-seq(model_sol$vec_date[2],model_sol$tstep,length=H)
  mylist<-list("P.tib"=P.tib,"r.tib"=r.tib,"date"=vec_date)
  return(mylist)
}

##Proposition 10: payoff on date t+h = exp(t(omega)%*%X)
#*H is the horizon of the bond [t+1:t+99], max H=98
#*omega.varphi =dim(X)*1 matrix
#*return r in percent

varphi<-function(model_sol,omega.varphi,H,X=model_sol$X,t=0){
  if((t+H)>=(model_sol$Tmax-1)){
    P.pi    <-model_sol$pi
    inf     <-rep(list(model_sol$inf_matx$pi.inf),t+H-(model_sol$Tmax-1)+1)
    P.pi    <-c(P.pi,inf)
    
    P.eta_1 <-model_sol$eta1
    inf     <-rep(list(model_sol$inf_matx$eta1.inf),t+H-(model_sol$Tmax-1)+1)
    P.eta_1 <-c(P.eta_1,inf)
    
    P.eta_0 <-model_sol$eta0
    inf     <-matrix(model_sol$inf_matx$eta0.inf,t+H-(model_sol$Tmax-1)+1,1)
    P.eta_0 <-rbind(P.eta_0,inf)
  }else{
    P.pi    <-model_sol$pi
    P.eta_1 <-model_sol$eta1
    P.eta_0 <-model_sol$eta0
  }
  
  #List of all our U.sh
  U.tk<-list(P.pi[[t+1]]+omega.varphi)
  P.a.pi<-matrix(NaN,H,1)
  for(h in 1:H){
    if((t+h)<=length(model_sol$pi)){
      P.a.pi[h]<-a1.fct.1(model_sol,P.pi[[t+h]],t+(h-1))
    }else{
      P.a.pi[h]<-a1.fct.inf(model_sol,P.pi[[t+h]]) 
    }
  }
  
  if(H>1){
    for (h in 2:H){
      uk<-matrix(NaN,nrow=model_sol$n.X,ncol=h)
      for(k in 1:(h-1)){
        if((t+k)<(model_sol$Tmax-1)){
          uk[,k] <--P.eta_1[[t+k+1]]-
            b1.fct.1(model_sol,P.pi[[t+k+1]],t+k)+
            P.pi[[t+k]]
        }else{
          uk[,k] <--P.eta_1[[t+k+1]]-
            b1.fct.inf(model_sol,P.pi[[t+k+1]])+
            P.pi[[t+k]]
        }
      }
      uk[,h]   <- P.pi[[t+h]]+omega.varphi
      U.tk[[h]]<-uk
    }
  }
  
  P.psi<-lapply(1:H,function(h){
    if(h==1){
      multi.lt.fct(model_sol,U.tk[[h]],h,X,t)
    }else{multi.lt.fct.Uh(model_sol,U.tk[[h]],X,t)}
    })
  
  varphi0<-matrix(NaN,H,1)
  varphi1<-list()
  
  for (h in 1:H){
    varphi0[h]  <--sum(P.eta_0[t+(1:h)])-
                   sum(P.a.pi[1:h])+
                   P.psi[[h]][[1]]
    if(t>(model_sol$Tmax-2)){
      varphi1[[h]]<--P.eta_1[[t+1]]-
        b1.fct.inf(model_sol,P.pi[[t+1]])+
        P.psi[[h]][[2]]
    }else{
      varphi1[[h]]<--P.eta_1[[t+1]]-
        b1.fct.1(model_sol,P.pi[[t+1]],t)+
        P.psi[[h]][[2]]
    }
  }
  
  varphi<-matrix(unlist(lapply(1:H,function(h)exp(varphi0[h]+
                                                    t(varphi1[[h]])%*%X))),
                 H,1)
  r.t<--log(varphi)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("varphi0"=varphi0,"varphi1"=varphi1,"P.t"=varphi,"r.t"=r.t,
               "date"=vec_date)
  return(mylist)
}

#Proposition 11: payoff on date t+h = exp(t(omega)%*%X)*1_{t(a)%*%X<b}
#*a is a vector dim(X)*1 and b is a scalar associated with the options 
#*a= payment associated with some components of X, and b=strike

varphi.hat<-function(model_sol,omega.v.hat,H,x,a,b,X=model_sol$X,t=0){
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  fx<-matrix(NaN,H,length(x))

  for (i in 1:length(x)){
    if(i==length(x)*round(i/length(x),1)){
      print(paste("Progress: ",100*i/length(x),"%",sep=""))
    }

    fx[,i]<-Im(varphi(model_sol,1i*a*x[i]+omega.v.hat,H,X,t)[[3]]*
                 exp(-1i*b*x[i]))/x[i]*dx[i]
  }
  print("done")
  varphi.hat<-varphi(model_sol,omega.v.hat,H,X,t)[[3]]/2-
              1/pi*matrix(apply(fx,1,sum),H,1)#X,t?

  r.hat<--log(varphi.hat)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.hat"=varphi.hat,"r.hat"=r.hat,"date"=vec_date)
  
  return(mylist)
}

#Corollary 3: payoff on date t+h = t(omega)%*%X

varphi.tilde<-function(model_sol,omega.v.tilde,H,X=model_sol$X,t=0){
  eps  <-10^-5
  o.ZCB<-matrix(0,model_sol$n.X,1)
  
  varphi.tilde<-(varphi(model_sol,eps*omega.v.tilde,H,X,t)[[3]]-
                   varphi(model_sol,o.ZCB,H,X,t)[[3]])/eps
  
  r.tilde<--log(varphi.tilde)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.tilde"=varphi.tilde,"r.tilde"=r.tilde,"date"=vec_date)
  
  return(mylist)
}

#Corollary 5: payoff on date t+h = t(omega)%*%X*1_{t(a)%*%X<b}
#*x is a grid that allow the function to do the approximated Fourier transform

varphi.bar<-function(model_sol,omega.v.bar,H,x,a,b,X=model_sol$X,t=0){
  eps       <-10^-5
  varphi.bar<-(varphi.hat(model_sol,eps*omega.v.bar,H,x,a,b,X,t)[[1]]-
                 varphi.hat(model_sol,0,H,x,a,b,X,t)[[1]])/eps
  
  print("done")
  r.bar<--log(varphi.bar)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.bar"=varphi.bar,"r.bar"=r.bar,"date"=vec_date)
  
  return(mylist)
}

################################################################################*FOURIER TRANSFORM
#Fourier set for a cdf for now (u=0) - CHECK FOR U!=0
#*i is the variable in X we are interested in
#*h is the maturity/horizon
#*u is set to 0 to compute the cdf
#*x is a grid for the integral
#*gamma is the set of possible values for the variable i

psi<-function(model_sol,u,h,i,X=model_sol$X){
  U    <-matrix(0,model_sol$n.X,length(u))
  U[i,]<-u
  psi  <-multi.lt.fct.N(model_sol,U,h,X)
  return(psi)
}

fourier<-function(model_sol,x,gamma,h,i,X=model_sol$X,u=0){
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  s1<-psi(model_sol,u+1i*x,h,i,X)
  fx<-outer(x,gamma,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
      dx[,1]
  f <-1/2-1/pi*apply(fx,2,sum)
  return(f)
}

################################################################################*MULTI LT, SIMPLE + COMPLETE
##Proposition 3: Multi-horizon Laplace transform of X
#*U is a matrix of dimension "dim(X)*H", with H>1, the maturity
#*return psi0, psi1 and psi_t.H

multi.lt.fct.Uh<-function(model_sol,Uh,X=model_sol$X,t=0){
  H  <-dim(Uh)[2]
  U  <-list(Uh[,H,drop=F])
  if((t+H)>(model_sol$Tmax-1)){
    a.h<-matrix(c(a1.fct.inf(model_sol,U[[1]]),rep(NaN,H-1)),H,1) 
    i<-1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i)>(model_sol$Tmax-2)&i<H){
      U[[i+1]]<-Uh[,H-i,drop=F]+b1.fct.inf(model_sol,U[[i]])
      a.h[i+1]<-a1.fct.inf(model_sol,U[[i+1]])
      i<-i+1
    }
    if(i<=(H-1)){
      for (k in i:(H-1)){
        U[[k+1]]<-Uh[,H-k,drop=F]+b1.fct.1(model_sol,U[[k]],t+H-k)
        a.h[k+1]<-a1.fct.1(model_sol,U[[k+1]],t+H-(k+1))
      } 
    }
  }else{
    a.h<-matrix(c(a1.fct.1(model_sol,U[[1]],t+H-1),rep(NaN,H-1)),H,1)   
    for (k in 1:(H-1)){
      U[[k+1]]<-Uh[,H-k,drop=F]+b1.fct.1(model_sol,U[[k]],t+H-k)
      a.h[k+1]<-a1.fct.1(model_sol,U[[k+1]],t+H-(k+1))
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct.1(model_sol,U[[H]],t))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct.1(model_sol,U[[H]],t) 
  }
  mylist<-list("psi.0"=psi.0,"psi.1"=psi.1,"uX_t.h"=uX_t.h)
  return(mylist)
}

##Corollary 1: Simple Multihorizon Laplace Transform
#*h corresponds to the maturity, max 99, then time-independent functions
#*U is a column vector of dimension "dim(X)*1"
#*return psi0, psi1 and psi_t.h

multi.lt.fct<-function(model_sol,U,h,X=model_sol$X,t=0){
  param<-model_sol$parameters
  
  U.h  <-list(U)
  a.sum<-0
  #recursive method with h going to time-indep maturity
  if((t+h)>(model_sol$Tmax-1)){
    i<-1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum     <-a.sum+a1.fct.inf(model_sol,U.h[[i]])
      U.h[[i+1]]<-b1.fct.inf(model_sol,U.h[[i]])
      
      i<-i+1
    }
    if(i<=h){
      for (k in i:h){
        a.sum     <-a.sum+a1.fct.1(model_sol,U.h[[k]],t+h-k)
        U.h[[k+1]]<-b1.fct.1(model_sol,U.h[[k]],t+h-k)
      }  
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }else{#Or max of time-dep maturity
    for (k in 1:h){
      a.sum     <-a.sum+a1.fct.1(model_sol,U.h[[k]],t+h-k)
      U.h[[k+1]]<-b1.fct.1(model_sol,U.h[[k]],t+h-k)
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  psi0<-a.sum
  psi1<-U.h[[h+1]]
  
  mylist<-list("psi0"=psi0,"psi1"=psi1,"uX_t.h"=uX_t.h)
  return(mylist)
}

##Corollary 1.N: Simple Multihorizon Laplace Transform
#*U is a matrix of dimension "dim(X)*N"
#*return psi_t.h

multi.lt.fct.N<-function(model_sol,U,h,X=model_sol$X,t=0){
  param    <-model_sol$parameters
  
  U.h  <-list(U)
  a.sum<-matrix(0,nrow=dim(U)[2])
  if((t+h)>(model_sol$Tmax-1)){
    i<-1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum     <-a.sum+a1.fct.inf(model_sol,U.h[[i]])
      U.h[[i+1]]<-b1.fct.inf(model_sol,U.h[[i]])
      
      i<-i+1
    }
    if(i<=h){
      for (k in i:h){
        a.sum     <-a.sum+a1.fct.N(model_sol,U.h[[k]],t+h-k)
        U.h[[k+1]]<-b1.fct.N(model_sol,U.h[[k]],t+h-k)
      } 
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
    
  }else{
    for (k in 1:h){
      a.sum     <-a.sum+a1.fct.N(model_sol,U.h[[k]],t+h-k)
      U.h[[k+1]]<-b1.fct.N(model_sol,U.h[[k]],t+h-k)
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  return(uX_t.h)
}

################################################################################*SDF SOLVING
##Proposition 8: Function that returns mu_u1 and mu_u0 @0 
#*General case with Tmax=100, end up with t=2015
#*If want a specific date--> mu_u.t.fct.all
#*return mu_u0.t, mu_u1.t

mu_u.t.fct<-function(model_sol){
  param    <-model_sol$parameters
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  
  mu_c1  <-model_sol$mu_c1
  mu_c0  <-model_sol$mu_c0
  mu_u1.t<-mu_u1.inf                                                            #mu_u1.99
  mu_u0.t<-mu_u0.inf                                                            #mu_u0.99
  #mu_u0 and mu_u1.98-0,99=inf.
  for (i in 2:model_sol$Tmax){
    b.t    <-b1.fct.1(model_sol,(1-param$gamma)*(mu_u1.t+mu_c1),model_sol$Tmax-i)
    a.t    <-a1.fct.1(model_sol,(1-param$gamma)*(mu_u1.t+mu_c1),model_sol$Tmax-i)
    mu_u1.t<-param$delta/(1-param$gamma)*b.t
    mu_u0.t<-param$delta*(mu_u0.t+mu_c0)+
             param$delta/(1-param$gamma)*a.t
  }
  
  mylist<-list("mu_u0.t"=mu_u0.t,"mu_u1.t"=mu_u1.t)
  return(mylist)
}

#*list with all mu_u1.t and mu_u0.t
#*from t=1 to t=98, inf(t=99) not in list
mu_u.t.fct.all<-function(model_sol){
  param    <-model_sol$parameters
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  Tmax     <-model_sol$Tmax

  mu_c1  <-model_sol$mu_c1
  mu_c0  <-model_sol$mu_c0
  
  mu_u1.t        <-list()
  mu_u1.t[[Tmax]]<-mu_u1.inf
  mu_u0.t        <-matrix(NaN,Tmax,1)
  mu_u0.t[Tmax]  <-mu_u0.inf
  for (i in 1:(Tmax-1)){
    mu_u1.t[[Tmax-i]]<-param$delta/(1-param$gamma)*
                         b1.fct.1(model_sol,
                                  (1-param$gamma)*(mu_u1.t[[Tmax-i+1]]+mu_c1),
                                   Tmax-i-1)
    mu_u0.t[Tmax-i]  <-param$delta*(mu_u0.t[Tmax-i+1]+mu_c0)+
                         param$delta/(1-param$gamma)*
                         a1.fct.1(model_sol,
                                  (1-param$gamma)*(mu_u1.t[[Tmax-i+1]]+mu_c1),
                                   Tmax-i-1)
  }
  #no inf
  mu_u1.t<-mu_u1.t[-Tmax]
  mu_u0.t<-mu_u0.t[-Tmax]
  #no 2015
  mu_u1  <-mu_u1.t[-1]
  mu_u0  <-mu_u0.t[-1]
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=Tmax-2)         #infinity not counted
  
  mylist<-list("mu_u0.t1"=mu_u0,"mu_u1.t1"=mu_u1,"date"=vec_date)
  return(mylist)
}

################################################################################*SOLVE MODEL, MEAN AND SIMUL
##Function solving the model for infinite mu=1
#*theta is a set of ini cond. for the optim. of mitig rate mu

model_solve<-function(model,theta,
                      indic_mitig=TRUE,mu.chosen=rep(param$mu0,model$Tmax)){
  #Data preparation
  model_sol<-model
  param    <-model$parameters
  tstep    <-model$tstep
  GN       <-0
  Tmax     <-model$Tmax
  #Determine the matrix of shocks                                               #+CHANGE in mu_dep
  omega.star.inf              <- matrix(0,model_sol$n.Z,model_sol$n.W)
  omega.star.inf[1,1]         <- param$sigma_a/(param$A_bar+1-param$delta_K)
  omega.star.inf[2,1]         <- param$sigma_a/(param$A_bar+1-param$delta_K)
  omega.star.inf[5,2]         <-(1-param$Phi[2,2]^2)^(1/2)*param$sigma_eta_f
  omega.star.inf[1,(n.eta+1)] <-0 #-1
  omega.star.inf[11,(n.eta+1)]<--1
  omega.star.inf[3,(n.eta+2)] <- 1
  omega.star.inf[14,3]        <- param$sigma_H
  if(indic_stocks==1){
    omega.star.inf[15,(n.eta+3)]<--1 
    omega.star.inf[16,]         <-t(model_sol$mu_pd.a0[[Tmax]]
                                    [(model_sol$n.Z+2):(model_sol$n.Z+model_sol$n.W+1)])
    omega.star.inf[19,]         <-t(model_sol$mu_r1.a1[[Tmax]]
                                    [(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)])
  }
  
  model_sol[["omega.star.inf"]]<-omega.star.inf
  
  ##############################################################################
  #Long term matrices and Gauss Newton Algorithm
  ##############################################################################
                                                                                #CHANGE vectorial rep.inf.
  A0.star.inf       <- diag(model_sol$n.Z)
  A0.star.inf[3,4]  <--1
  A0.star.inf[5,6]  <--param$tau/(log(2)*param$m_pi)*1/(param$m0)
  A0.star.inf[9,5]  <--param$xi_1
  A0.star.inf[12,3] <--1
  A0.star.inf[13,1] <--1
  A0.star.inf[14,9] <--model_sol$tstep*param$a_sat-param$b_sat
  A0.star.inf[1,14] <- 0 #param$b_sk
  if(indic_stocks==1){
    A0.star.inf[15,1] <--param$chi.div
    A0.star.inf[16,]  <--t(model_sol$mu_pd.a0[[Tmax]][2:(model_sol$n.Z+1)])
    A0.star.inf[16,16]<-1+A0.star.inf[16,16]
    A0.star.inf[17,15]<--1
    A0.star.inf[17,16]<--model_sol$kap1.a0[Tmax]
    A0.star.inf[18,17]<--1
    
  }

  varphi           <- matrix(0,model_sol$n.eta,model_sol$n.eta)
  varphi[1,1]      <- param$varphi_11
  varphi[2,1]      <- param$varphi_12
  varphi[1,2]      <- param$varphi_21
  varphi[2,2]      <- param$varphi_22
  varphi[3,2]      <- param$varphi_23
  varphi[2,3]      <- param$varphi_32
  varphi[3,3]      <- param$varphi_33
  
  A1.star.inf           <- matrix(0,model_sol$n.Z,model_sol$n.Z)
  A1.star.inf[2,2]      <- 1
  A1.star.inf[4,2]      <- 0
  A1.star.inf[6,3]      <- model_sol$tstep/3.666
  A1.star.inf[6:8,6:8]  <- varphi%^%(model_sol$tstep)                           #CHANGE if add eta
  A1.star.inf[9,9]      <- 1-param$xi_1*(param$lambda+param$xi_2)
  A1.star.inf[9,10]     <- param$xi_1*param$xi_2
  A1.star.inf[10,9]     <- param$xi_3
  A1.star.inf[10,10]    <- 1-param$xi_3
  A1.star.inf[11,11]    <- 1
  A1.star.inf[12,12]    <- 1
  A1.star.inf[13,13]    <- 1
  A1.star.inf[14,14]    <- 1
  A1.star.inf[14,9]     <--param$b_sat
  A1.star.inf[1,14]     <- param$b_sk
  if(indic_stocks==1){
    A1.star.inf[15,15]    <- 1
    A1.star.inf[17,15]    <--1
    A1.star.inf[17,16]    <--1
    A1.star.inf[18,18]    <- 1
    A1.star.inf[19,]      <- t(model_sol$mu_r1.a1[[Tmax]][1:model_sol$n.Z])
  }

  
  omega0.star.inf      <- matrix(0,model_sol$n.Z,1)
  omega0.star.inf[1,1] <- log(param$delta)+log(param$A_bar+1-param$delta_K)
  omega0.star.inf[3,1] <- 0
  omega0.star.inf[4,1] <- 0
  omega0.star.inf[5,1] <- param$tau/log(2)*(log(param$m0)-1)+param$phi_1
  omega0.star.inf[14,1]<--model_sol$tstep*param$a_sat*param$T_0s+param$c_sat
  if(indic_stocks==1){
    omega0.star.inf[15,1]<-model_sol$mu_div.0[Tmax]
    omega0.star.inf[16,1]<-model_sol$mu_pd.a0[[Tmax]][1]
    omega0.star.inf[17,1]<-model_sol$kap0.om0[Tmax]
    omega0.star.inf[19,1]<-model_sol$mu_r0.om0[Tmax]
  }

  
  A1.inf      <-solve(A0.star.inf)%*%A1.star.inf
  omega0.inf  <-solve(A0.star.inf)%*%omega0.star.inf
  omega.inf   <-solve(A0.star.inf)%*%omega.star.inf
  
  model_sol[["varphi"]]       <-varphi
  model_sol[["A0.star.inf"]]  <-A0.star.inf
  model_sol[["A1.inf"]]       <-A1.inf
  model_sol[["omega0.inf"]]   <-omega0.inf
  model_sol[["omega.inf"]]    <-omega.inf

  ##Infinite components (t+100)
  # solution initial values
  betawGN<- cbind(matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z),          #CHANGE IF ADD GAMMA0
               rbind(
                 matrix(0,model_sol$n.Z,model_sol$n.eta),
                 t(model_sol$parameters$Phi),
                 matrix(0,model_sol$n.W-model_sol$n.eta,model_sol$n.eta)
               ),
               model_sol$parameters$mu_d*model_sol$parameters$ell1.D,
               model_sol$parameters$mu_n*model_sol$ell1.N.tilde,
               if(indic_stocks==1){
                 model_sol$parameters$mu_div*model_sol$param$ell1.div 
               }
               )
  betaGN <- cbind(rbind(t(model_sol$A1.inf),
                        matrix(0,model_sol$n.W,model_sol$n.Z)),
                  matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.W)
                  )+
    betawGN%*%rbind(matrix(0,model_sol$n.Z,model_sol$n.Z+model_sol$n.W),
                    cbind(t(model_sol$omega.inf),diag(model_sol$n.W))
                    )
  model_sol[["M"]]<-betaGN
  muc1GN  <- rbind(1,matrix(0,model_sol$n.Z+model_sol$n.W-1))
  muu1_ini<- model_sol$parameters$delta*
    solve(diag(model_sol$n.Z+model_sol$n.W)-model_sol$parameters$delta*betaGN)%*%
    betaGN%*%muc1GN
  GN      <-muu1_ini

  Newton<-GN.function(model_sol,GN)
  
  mu_c1 <-Newton$mu_c1
  mu_c0 <-0
  mu_u1 <-Newton$mu_u1
  mu_u0 <-param$delta/(1-param$delta)*
          mu_c0+param$delta/((1-param$delta)*(1-param$gamma))*
          a1.fct.inf(model_sol,(1-param$gamma)*(mu_u1+mu_c1))
  
  model_sol[["mu_c1"]]<-mu_c1
  model_sol[["mu_c0"]]<-mu_c0
  model_sol[["mu_u1"]]<-mu_u1
  model_sol[["mu_u0"]]<-mu_u0
  model_sol[["ite"]]  <-Newton$ite
  model_sol[["dev"]]  <-Newton$dev
  
  ##############################################################################
  #Optimization of mu - X and exogenous vectors
  ##############################################################################
  #####
  #Exogenous variables
  #####
  f_ex  <-matrix(NaN, nrow=Tmax, 1)
  E_land<-matrix(NaN, nrow=Tmax, 1)
  gsigma<-matrix(NaN, nrow=Tmax, 1)
  sigma <-matrix(NaN, nrow=Tmax, 1)
  mu    <-matrix(NaN, nrow=Tmax, 1)
  bp    <-matrix(NaN, nrow=Tmax, 1)
  bc    <-matrix(NaN, nrow=Tmax, 1)
  
  #Radiative forcings
  f_ex               <- matrix(rep(param$phi_0,Tmax),Tmax,1)                        
  f_ex[1:17]         <- f_ex[1:17]+(1/17)*(param$phi_1-param$phi_0)*((1:17)-1)
  f_ex[18:Tmax]      <- f_ex[18:Tmax] + (param$phi_1-param$phi_0)
  
  #Emissions from deforestation
  E_land[1:Tmax]<- param$eps_0*(1-param$rho)**((1:(Tmax))-1)           
  
  #Carbon intensity
  gsigma[1]<-param$gsigma1                                                      
  for(i in 2:Tmax) gsigma[i] <- gsigma[i-1]*((1+param$delsigma)**tstep)    
  sigma[1] <- param$sigma0                                                         
  for(i in 2:Tmax) sigma[i]  <- sigma[i-1] * exp(gsigma[i-1] * tstep)      
  
  #Abatement cost exogenous components
  bp[1:Tmax]<-param$pback*(1-param$gback)**((1:Tmax)-1)              
  bc[1:Tmax]<-bp[1:Tmax]*sigma[1:Tmax]/param$theta2/1000
  
  #####
  #Initial states in 2015                                                       #CHANGE STATES INI COND.
  #####
  Z    <- matrix(NaN,model_sol$n.Z,1)
  Z[1] <- log(param$delta)+log((1-bc[1]*param$mu0**(param$theta2))*param$A_bar+
                                 1-param$delta_K)
  Z[2] <- param$ini_tildey
  Z[3] <- param$ini_E
  Z[4] <- param$ini_Eind
  Z[5] <-(param$ini_Mat/param$m_pi-1)*param$tau/log(2)+f_ex[1]
  Z[6] <- param$ini_Mat
  Z[7] <- param$ini_Mup
  Z[8] <- param$ini_Mlo
  Z[9] <- param$ini_Tat 
  Z[10]<- param$ini_Tlo
  Z[11]<- param$ini_CumD
  Z[12]<- Z[3]
  Z[13]<- 0
  Z[14]<- param$ini_H
  if(indic_stocks==1){
    Z[15]<- 0
    Z[16]<- model_sol$ini_pd
    Z[17]<- 0
    Z[18]<- 0
    Z[19]<- 0
  }
  W    <- matrix(0,model_sol$n.W,1)
  #t=0, value of the state variables at time 0
  #X=[delc,tilde_y,E,E_ind,Forc,M_at,M_up,M_lo,T_at,T_lo,Cum_D,Cum_E,Cum_dc,W]
  X<-rbind(Z,W)
  
  model_sol$parameters$ini_delc   <-Z[1]
  model_sol$parameters$ini_Cumdelc<-Z[1]
  
  #Keep all elements st [1] is t for t=0
  model_sol[["f_ex"]]  <-f_ex
  model_sol[["E_land"]]<-E_land
  model_sol[["gsigma"]]<-gsigma
  model_sol[["sigma"]] <-sigma
  model_sol[["bp"]]    <-bp
  model_sol[["bc"]]    <-bc
  
  model_sol[["X"]]     <-X
  model_sol[["n.X"]]   <-length(X)
  
  #####
  #Optimization calculations
  #####
  if (indic_mitig){
    opt<-list()
    for (i in 1:length(theta)){
      opt[[i]] <- res.optim(model_sol,theta[[i]])
    }
    maxi     <-min(unlist(extract(opt,2)))
    best     <-which(maxi==unlist(extract(opt,2)))[1]
    theta.opt<-unlist(opt[[best]][1])
    mu       <-mu.function(model_sol,theta.opt)
    
    
    model_sol[["theta.opt"]]<-theta.opt
    model_sol[["u0"]]       <-abs(maxi)
  }else{
    mu   <- mu.chosen
  }
  
  model_sol[["mu"]]   <-mu
  #Date 0 storage
  list2015             <-list()
  
  ##############################################################################
  #Exogenous variables depending on mu optimized (or chosen)
  #####
  model_matrix<-mu_dep(model_sol,model_sol$mu)
  #####
 
  #Keep all elements st [1] is t for t=0
  model_sol[["AC"]]      <-model_matrix$AC
  model_sol[["kappa"]]   <-model_matrix$kappa
  
  #remove first period st starting from t=0, first element is t=1
  model_sol[["A1"]]      <-model_matrix$A1[-1]
  model_sol[["omega"]]   <-model_matrix$omega[-1]
  model_sol[["omega0"]]  <-model_matrix$omega0[-1]
  
  #Date 0 storage
  list2015[["A1"]]       <-model_matrix$A1[[1]]
  list2015[["omega"]]    <-model_matrix$omega[[1]]
  list2015[["omega0"]]   <-model_matrix$omega0[[1]]

  if (indic_mitig==FALSE){
    mu_u1.t  <-model_sol$mu_u1
    mu_u0.t  <-model_sol$mu_u0
    mu_c1    <-model_sol$mu_c1
    mu_c0    <-model_sol$mu_c0

    mu_u.1<-mu_u.t.fct(model_sol)
    
    if(is.na(mu_u.1[[1]])){
      u0 <- 10000
      model_sol[["u0"]] <-u0
    }else{
      u0<-log(model_sol$c0)+mu_u.1[[1]]+t(mu_u.1[[2]])%*%model_sol$X
      model_sol[["u0"]] <-u0
    }
  }


  #####
  #*construct mu_u1.t
  #####
  #Proposition 9: SDF construction
  P.mu_u.t1 <-mu_u.t.fct.all(model_sol)
  P.mu_u1.t1<-P.mu_u.t1[[2]]

  P.pi  <-lapply(1:(Tmax-2),                                                    
                 function(t){(1-param$gamma)*P.mu_u1.t1[[t]]-
                     param$gamma*mu_c1})
  
  #####
  #*construct eta.0                                                             
  #####
  P.eta_0  <-matrix(unlist(lapply(1:(length(P.pi)),function(t){
              -log(param$delta)+
                a1.fct.1(model_sol,(1-param$gamma)*(P.mu_u1.t1[[t]]+mu_c1),t-1)-
                a1.fct.1(model_sol,P.pi[[t]],t-1)})),                           #t=0 for a.t
              length(P.pi),1)
  
  #####
  #*construct eta.1                                                             
  #####
  P.eta_1  <-lapply(1:(length(P.pi)),
                    function(t){b1.fct.1(model_sol,
                                         (1-param$gamma)*(P.mu_u1.t1[[t]]+mu_c1),
                                         t-1)-
                                b1.fct.1(model_sol,P.pi[[t]],t-1)})             #t=0 for b.t
  
  
  model_sol[["mu_u1.t1"]]<-P.mu_u1.t1
  model_sol[["pi"]]      <-P.pi
  model_sol[["eta1"]]    <-P.eta_1
  model_sol[["eta0"]]    <-P.eta_0
  
  #Infinite economy:
  inf_matx<-list("mu_u1"=model_sol$mu_u1)
  
  inf_matx[["pi.inf"]]  <-(1-param$gamma)*inf_matx$mu_u1-param$gamma*mu_c1
  inf_matx[["eta0.inf"]]<--log(param$delta)+
    a1.fct.inf(model_sol,(1-param$gamma)*(inf_matx$mu_u1+mu_c1))-
    a1.fct.inf(model_sol,inf_matx$pi.inf)
  inf_matx[["eta1.inf"]]<-b1.fct.inf(model_sol,
                                     (1-param$gamma)*(inf_matx$mu_u1+mu_c1))-
    b1.fct.inf(model_sol,inf_matx$pi.inf)
  

  
  model_sol$ini_matx<-c(list2015,model_sol$ini_matx)
  model_sol$inf_matx<-c(inf_matx,model_sol$inf_matx)
  ##############################################################################
  return(model_sol)
}
#####
#Optimization of mu - Functions
#####
#Function to prevent agents to mitigate too quickly
mu.function <- function(model_sol,theta,
                        date.min.mu.equal.1=passive.mu,max.mu.date.0=.2){
  a <- -log(max.mu.date.0) + abs(theta[1])
  b <- abs(theta[2])
  
  if(a<b*date.min.mu.equal.1){
    b <- a/date.min.mu.equal.1
  }
  mu  <-pmin(exp(-a+b*(1:model_sol$Tmax)),1)
  return(mu)
}


#Compute the vectorial repr. of the model + exogenous dependent fct(mu)
mu_dep<-function(model_sol,mu){
  param     <-model_sol$parameters
  Tmax      <-model_sol$Tmax

  AC    <-matrix(NaN, nrow=Tmax, 1)
  kappa <-matrix(NaN, nrow=Tmax, 1)
  
  #Exogenous equations independent of mu
  #rad. forcings
  f_ex   <- model_sol$f_ex                        
  #emissions from deforestation
  E_land <- model_sol$E_land           
  #Carbon intensity
  gsigma <- model_sol$gsigma
  sigma  <- model_sol$sigma
  #Abatement cost exogenous components
  bp     <- model_sol$bp            
  bc     <- model_sol$bc
  
  #Exogenous equations dependent on mu
  #Abatement Cost
  AC[1]        <-bc[1]*param$mu0**(param$theta2)                                   
  AC[2:Tmax]   <-bc[2:Tmax]*mu[2:Tmax]**param$theta2  

  #kappa
  kappa[1]      <- sigma[1]*(1-param$mu0)*param$q0*
    exp((log(param$delta)+log((1-AC[1])*param$A_bar+1-param$delta_K)+
             1/2*((1-AC[1])*param$sigma_a/
                    ((1-AC[1])*param$A_bar+1-param$delta_K))^2))
  
  for (i in 2:Tmax){
    kappa[i] <-(1-mu[i])*sigma[i]*param$q0*
      exp(sum(log(param$delta)+log((1-AC[1:i])*param$A_bar+1-param$delta_K)+
                1/2*((1-AC[1:i])*param$sigma_a/
                       ((1-AC[1:i])*param$A_bar+1-param$delta_K))^2))
  }
  
  
  #Vectorial representation                                                     #CHANGE MATX.(non inf)
  if(indic_stocks==1){
    A0.star<-list()
    for (i in 1:Tmax){
      A0_i        <- diag(model_sol$n.Z)
      A0_i[3,4]   <--1
      A0_i[5,6]   <--param$tau/(log(2)*param$m_pi)*1/(param$m0)
      A0_i[9,5]   <--param$xi_1
      A0_i[12,3]  <--1
      A0_i[13,1]  <--1
      A0_i[14,9]  <--model_sol$tstep*param$a_sat-param$b_sat
      A0_i[1,14]  <- param$b_sk
      A0_i[15,1]  <- -param$chi.div
      A0_i[16,]   <--t(model_sol$mu_pd.a0[[i]][2:(model_sol$n.Z+1)])
      A0_i[16,16] <- 1+A0_i[16,16]
      A0_i[17,15] <--1 
      A0_i[17,16] <--model_sol$kap1.a0[i]
      A0_i[18,17] <--1 
      A0.star[[i]]<- A0_i
    } 
  }else{
    A0.star<-model_sol$A0.star.inf
  }
  
  A1.star<-list()
  for (i in 1:Tmax){
    A1_i           <- matrix(0,model_sol$n.Z,model_sol$n.Z)
    A1_i[2,2]      <- 1
    A1_i[4,2]      <- kappa[i]
    A1_i[6,3]      <- model_sol$tstep/3.666
    A1_i[6:8,6:8]  <- model_sol$varphi%^%(model_sol$tstep)
    A1_i[9,9]      <- 1-param$xi_1*(param$lambda+param$xi_2)
    A1_i[9,10]     <- param$xi_1*param$xi_2
    A1_i[10,9]     <- param$xi_3
    A1_i[10,10]    <- 1-param$xi_3
    A1_i[11,11]    <- 1
    A1_i[12,12]    <- 1
    A1_i[13,13]    <- 1
    A1_i[14,14]    <- 1
    A1_i[14,9]     <--param$b_sat
    A1_i[1,14]     <- param$b_sk
    if(indic_stocks==1){
      A1_i[15,15]    <- 1
      A1_i[17,15]    <--1
      A1_i[17,16]    <--1
      A1_i[18,18]    <- 1
      A1_i[19,]      <- model_sol$mu_r1.a1[[i]][1:model_sol$n.Z]
    }
    A1.star[[i]]   <- A1_i 
    
  }
  
  omega0.star<-list()
  for(i in 1:Tmax){
    omega0_i        <- matrix(0,model_sol$n.Z,1)
    if(indic_tp==1){
      omega0_i[1,1]   <- model_sol$mu_c[i]
    }else{
      omega0_i[1,1]   <- log(param$delta)+
                          log((1-AC[i])*param$A_bar+1-param$delta_K)
    }
    omega0_i[3,1]   <- E_land[i]
    omega0_i[4,1]   <- kappa[i]
    omega0_i[5,1]   <- param$tau/log(2)*(log(param$m0)-1)+f_ex[i]
    omega0_i[14,1]  <--model_sol$tstep*param$a_sat*param$T_0s+param$c_sat
    if(indic_stocks==1){
      omega0_i[15,1]  <- model_sol$mu_div.0[i]
      omega0_i[16,1]  <- model_sol$mu_pd.a0[[i]][1]
      omega0_i[17,1]  <- model_sol$kap0.om0[i]
      omega0_i[19,1]  <- model_sol$mu_r0.om0[i]
    }
    omega0.star[[i]]<- omega0_i
  }
  
  omega.star<-list()
  for(i in 1:Tmax){
    omega_i              <- matrix(0,model_sol$n.Z,model_sol$n.W)
    omega_i[1,1]         <- param$sigma_a*
                              (1-AC[i])/((1-AC[i])*param$A_bar+1-param$delta_K)
    omega_i[2,1]         <- param$sigma_a*
                              (1-AC[i])/((1-AC[i])*param$A_bar+1-param$delta_K)
    omega_i[5,2]         <-(1-param$Phi[2,2]^2)^(1/2)*param$sigma_eta_f
    omega_i[1,(n.eta+1)] <--1
    omega_i[11,(n.eta+1)]<--1
    omega_i[3,(n.eta+2)] <- 1
    omega_i[14,3]        <- param$sigma_H
    if(indic_stocks==1){
      omega_i[15,(n.eta+3)]<--1 
      omega_i[16,]         <-t(model_sol$mu_pd.a0[[i]]
                               [(model_sol$n.Z+2):(model_sol$n.Z+model_sol$n.W+1)])
      omega_i[19,]         <-t(model_sol$mu_r1.a1[[i]]
                               [(model_sol$n.Z+1):model_sol$n.X])
    }
    omega.star[[i]]      <- omega_i
  }
  if(indic_stocks==1){
    A1      <-lapply(1:Tmax,function(i)solve(A0.star[[i]])%*%
                       A1.star[[i]])
    omega0  <-lapply(1:Tmax,function(i)solve(A0.star[[i]])%*%
                       omega0.star[[i]])
    omega   <-lapply(1:Tmax,function(i)solve(A0.star[[i]])%*%
                       omega.star[[i]])
  }else{
    A1      <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                       A1.star[[i]])
    omega0  <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                       omega0.star[[i]])
    omega   <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                       omega.star[[i]]) 
  }
  
  mylist  <-list("A1"=A1,"omega0"=omega0,"omega"=omega,
                 "AC"=AC,"kappa"=kappa)
  return(mylist)
}

#######################################
#Results of optimization problem of mu, linked to utility function
#return the value of optimization and the a/b (theta1/2)
res.optim <-function(model_sol,theta){
  #print(utility.optim(model_sol,theta))
  res.optim<-optim(par=theta,utility.optim,
                   model_sol=model_sol,
                   gr = NULL,
                   method="Nelder-Mead",
                   #method="CG",
                   #method="BFGS",
                   control=list(trace=F,maxit=model_sol$MAXIT),
                   hessian=FALSE)
  mylist<-list("res.optim$par"=res.optim$par,
               "res.optim$value"=res.optim$value)
  return(mylist)
}
#construction of u0 for the optimization of mitigation rate mu
utility.optim<-function(model_sol,theta){
  param     <-model_sol$parameters
  c.t       <-model_sol$c0
  tstep     <-model_sol$tstep
  omega.star<-model_sol$omega.star
  Tmax      <-model_sol$Tmax
  
  #utility
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  mu_c1    <-model_sol$mu_c1
  mu_c0    <-model_sol$mu_c0                                                   
  mu_u1.t  <-mu_u1.inf                                                            
  mu_u0.t  <-mu_u0.inf
  
  #Parameters to optimize
  a <-theta[1]
  b <-theta[2]
  
  #List of variables
  mu<-matrix(NaN, nrow=Tmax, 1)
  
  #Emissions control rate
  mu[1:Tmax]<-mu.function(model_sol,c(a,b))
  
  model_matrix<-mu_dep(model_sol,mu)
  omega0      <-model_matrix$omega0
  omega       <-model_matrix$omega
  A1          <-model_matrix$A1
  
  model_sol[["A1"]]    <-A1[-1]
  model_sol[["omega"]] <-omega[-1]
  model_sol[["omega0"]]<-omega0[-1]
  
  mu_u.1<-mu_u.t.fct(model_sol) 
  #if NaN, returns high utility s.t. I don't mitigate and big loss
  if(is.na(mu_u.1[[1]])){
    return(-10000)
  }
  
  u0<-log(model_sol$c0)+mu_u.1[[1]]+t(mu_u.1[[2]])%*%model_sol$X
  return(-u0)
}

################################################################################
##Function to compute theoretical mean and variance
#*h is the end date of estimations, \in(1,99)
#*list starts in start_date+tstep

EV.fct<-function(model_sol,h=NaN){
  if(is.na(h)){
    t <- length(model_sol$vec_date)
  }else{
    t <- h
  }
  
  param <-model_sol$parameters
  if(t>(model_sol$Tmax-1)){
    omega0<-model_sol$omega0
    inf   <-rep(list(model_sol$omega0.inf),t-model_sol$Tmax+1)
    omega0<-c(omega0,inf)
    
    omega <-model_sol$omega
    inf   <-rep(list(model_sol$omega.inf),t-model_sol$Tmax+1)
    omega <-c(omega,inf)
    
    A1    <-model_sol$A1
    inf   <-rep(list(model_sol$A1.inf),t-model_sol$Tmax+1)
    A1    <-c(A1,inf)
  }else{
    omega0<-model_sol$omega0
    omega <-model_sol$omega
    A1    <-model_sol$A1
  }
  
  X     <-model_sol$X
  n.Z   <-model_sol$n.Z
  n.W   <-model_sol$n.W
  n.eta <-model_sol$n.eta
  
  #Proposition 4+5: Conditional mean and variance of W
  #alpha1.w, beta1.w, alpha2.w, beta2.w                                         #CHANGE if add gamma0
  alpha1.w                                <-matrix(0,model_sol$n.X,1)
  alpha1.w[(n.Z+1)+n.eta]                 <-param$mu_d*param$ell0.D
  alpha1.w[(n.Z+2)+n.eta]                 <-param$mu_n*param$ell0.N
  if(indic_stocks==1){
    alpha1.w[(n.Z+3)+n.eta]                 <-param$mu_div*param$ell0.div 
  }
  
  beta1.w                                 <-matrix(0,model_sol$n.X,model_sol$n.X)
  beta1.w[(n.Z+1):(n.Z+n.eta),
          (n.Z+1):(n.Z+n.eta)]            <-param$Phi                             
  beta1.w[(n.Z+1+n.eta),]                 <-param$mu_d*t(param$ell1.D)
  beta1.w[(n.Z+2+n.eta),]                 <-param$mu_n*t(model_sol$ell1.N.tilde)
  if(indic_stocks==1){
    beta1.w[(n.Z+3+n.eta),]                 <-param$mu_div*t(param$ell1.div) 
  }

  
  alpha2.w                                <-matrix(0,n.W,n.W)
  alpha2.w[1:n.eta,1:n.eta]               <-diag(n.eta)
  alpha2.w[(1+n.eta),(1+n.eta)]           <-2*param$mu_d^2*param$ell0.D
  alpha2.w[(2+n.eta),(2+n.eta)]           <-2*param$mu_n^2*param$ell0.N
  if(indic_stocks==1){
    alpha2.w[(3+n.eta),(3+n.eta)]           <-2*param$mu_div^2*param$ell0.div 
  }
  
  
  beta2.w.D                               <-matrix(0,n.W,n.W)
  beta2.w.D[(1+n.eta),(1+n.eta)]          <-2*param$mu_d^2
  beta2.w.N                               <-matrix(0,n.W,n.W)
  beta2.w.N[(2+n.eta),(2+n.eta)]          <-2*param$mu_n^2
  if(indic_stocks==1){
    beta2.w.N[(3+n.eta),(3+n.eta)]          <-2*param$mu_div^2 
  }
  
  alpha2.w<-matrix(alpha2.w,n.W*n.W,1)
  beta2.w <-matrix(beta2.w.D,n.W*n.W,1)%*%t(param$ell1.D)+
            matrix(beta2.w.N,n.W*n.W,1)%*%t(model_sol$ell1.N.tilde)
  remove(beta2.w.D,beta2.w.N)
  
  EV<-list()
  
  #Proposition 6+7: Conditional mean and variance of X
  #alpha1.k1
  alpha1.k1<-lapply(1:t,function(i){
    rbind(omega0[[i]],matrix(0,model_sol$n.W,1))+
      cbind(matrix(0,model_sol$n.X,model_sol$n.Z),
            rbind(omega[[i]],diag(model_sol$n.W)))%*%alpha1.w
  })
  
  #beta1.k1
  beta1.k1<-lapply(1:t,function(i){
    cbind(rbind(A1[[i]],matrix(0,model_sol$n.W,model_sol$n.Z)),
          matrix(0,model_sol$n.X,model_sol$n.W))+
      cbind(matrix(0,model_sol$n.X,model_sol$n.Z),
            rbind(omega[[i]],diag(model_sol$n.W)))%*%beta1.w
  })
  
  #Conditional mean for all h
  #alpha1.tk
  alpha1.tk<-list(alpha1.k1[[1]])
  if(t>1){
    for (i in 2:t){
      alpha1.tk[[i]]<-alpha1.k1[[i]]+beta1.k1[[i]]%*%(alpha1.tk[[i-1]])
    }
  }
  
  #beta1.tk
  beta1.tk<-list(beta1.k1[[1]])
  if(t>1){
    for (i in 2:t){
      beta1.tk[[i]]<-beta1.k1[[i]]%*%(beta1.tk[[i-1]])
    }
  }
  
  
  EV[["beta1.k1"]]<-beta1.k1
  
  #####
  #EXh#
  #####
  EXh<-list()
  for (i in 1:t){
    EXh[[i]]<-alpha1.tk[[i]]+beta1.tk[[i]]%*%X
  }
  
  EX<-lapply(1:model_sol$n.X,function(i){extract(EXh,i)})                       #CHANGE if add var
   
  names(EX)     <- c("delc",
                     "y_tilde",
                     "E",
                     "E_ind",
                     "Forc",
                     "M_at",
                     "M_up",
                     "M_lo",
                     "T_at",
                     "T_lo",
                     "Cum_D",
                     "Cum_E",
                     "Cum_dc",
                     "H",
                     if(indic_stocks==1){
                       c("Cum_div",
                         "pd",
                         "r.s",
                         "Cum_r.s",
                         "r.1")
                     },
                     "eta.A",
                     "eta.F",
                     "eta.H",
                     "D",
                     "N",
                     if(indic_stocks==1){
                       "D.div" 
                     }
  )
  
  EV[["EX"]]<-EX
  
  ######################
  #Conditional Variance#
  ######################
  
  #alpha2.k1
  alpha2.k1<-lapply(1:t,function(i){
    ((t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
        rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (omega[[i]]%x%omega[[i]])+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (diag(model_sol$n.W)%x%omega[[i]])+
       (t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       (omega[[i]]%x%diag(model_sol$n.W))+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       diag(model_sol$n.W*model_sol$n.W))%*%alpha2.w})
  
  #beta2.k1
  beta2.k1<-lapply(1:t,function(i){
    ((t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
        rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (omega[[i]]%x%omega[[i]])+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (diag(model_sol$n.W)%x%omega[[i]])+
       (t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       (omega[[i]]%x%diag(model_sol$n.W))+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       diag(model_sol$n.W*model_sol$n.W))%*%beta2.w})
  
  #####
  #VXh#
  #####
  
  #alpha2.tk
  alpha2.tk<-list(alpha2.k1[[1]])
  if(t>1){
    for (i in 2:t){
      alpha2.tk[[i]]<-alpha2.k1[[i]]+beta2.k1[[i]]%*%alpha1.tk[[i]]+
        (beta1.k1[[i]]%x%beta1.k1[[i]])%*%alpha2.tk[[i-1]]
    }
  }
  
  
  #beta2.tk
  beta2.tk<-list(beta2.k1[[1]])
  if(t>1){
    for (i in 2:t){
      beta2.tk[[i]]<-beta2.k1[[i]]%*%beta1.tk[[i]]+
        (beta1.k1[[i]]%x%beta1.k1[[i]])%*%beta2.tk[[i-1]]
    }
  }
 
  
  
  #vecVXh
  vecVXh  <-lapply(1:t,
                   function(x)alpha2.tk[[x]]+beta2.tk[[x]]%*%X)
  
  VXh     <-lapply(1:t,
                   function(x)matrix(vecVXh[[x]],
                                     model_sol$n.X,model_sol$n.X))
  
  VXh.diag<-lapply(1:t,
                   function(x)diag(VXh[[x]]))
  
  VX      <-lapply(1:model_sol$n.X, 
                   function(i)extract(VXh.diag,i))
  
  #lower and upper bound
  upper   <-lapply(1:model_sol$n.X,function(x){EX[[x]]+2*(VX[[x]])^(1/2)})
  lower   <-lapply(1:model_sol$n.X,function(x){EX[[x]]-2*(VX[[x]])^(1/2)})
  
  bounds  <-lapply(1:model_sol$n.X,function(x){rbind(lower[[x]],upper[[x]])})
  
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=t)
  
  EV[["CovX"]]  <-VXh
  EV[["VX"]]    <-VX
  EV[["bounds"]]<-bounds
  EV[["date"]]  <-vec_date
  return(EV)
}

################################################################################
                                                                               #CHANGE
##Function to simulate state variables of our model
#*New approximation of the radiative forcings (linear)
#*nb.simul is number of periods for 1 simulation \in(1,99)
#*nb.traj is the number of simulations done for nb.simul periods
#*for t>99, time-independent matrices

simul.function<-function(model_sol,nb.simul.t,nb.traj){
  #useful parameters
  nb.simul<-nb.simul.t+1 #take into account the t=0
    
  tstep<-model_sol$tstep
  t    <-1:nb.simul
  param<-model_sol$parameters
  
  if(nb.simul>model_sol$Tmax){
    omega0<-model_sol$omega0
    inf   <-rep(list(model_sol$omega0.inf),nb.simul-model_sol$Tmax)
    omega0<-c(omega0,inf)
    
    omega <-model_sol$omega
    inf   <-rep(list(model_sol$omega.inf),nb.simul-model_sol$Tmax)
    omega <-c(omega,inf)
    
    A1    <-model_sol$A1
    inf   <-rep(list(model_sol$A1.inf),nb.simul-model_sol$Tmax)
    A1    <-c(A1,inf)
  }else{
    omega0<-model_sol$omega0
    omega <-model_sol$omega
    A1    <-model_sol$A1
  }
  
  #Shocks
  eta  <-matrix(0, nrow=model_sol$n.eta, ncol=nb.traj)
  D    <-matrix(0, nrow=nb.simul, ncol=nb.traj)
  N    <-matrix(0, nrow=nb.simul, ncol=nb.traj)
  if(indic_stocks==1){
    D.div<-matrix(0, nrow=nb.simul, ncol=nb.traj) 
  }
  
  ##endogenous equations
  #Initial conditions
  delc   <-model_sol$X[1]                                                       
  y_tilde<-model_sol$X[2]
  E      <-model_sol$X[3]                                                       #CDICE2021
  E_ind  <-model_sol$X[4]                                                       #CDICE2021
  Forc   <-model_sol$X[5]                                                       #CDICE2021
  M_at   <-model_sol$X[6]                                                       #CDICE2021,mateq0
  M_up   <-model_sol$X[7]                                                       #CDICE2021,mueq0
  M_lo   <-model_sol$X[8]                                                       #CDICE2021,mleq0
  T_at   <-model_sol$X[9]                                                       #CDICE2021,tat0
  T_lo   <-model_sol$X[10]                                                      #CDICE2021,tlo0
  Cum_D  <-model_sol$X[11]
  Cum_E  <-model_sol$X[12]
  Cum_dc <-model_sol$X[13]
  H      <-model_sol$X[14]
  if(indic_stocks==1){
    Cum_div<-model_sol$X[15]
    pd     <-model_sol$X[16]
    r.s    <-model_sol$X[17]
    Cum_r.s<-model_sol$X[18]
    r.1    <-model_sol$X[19]
  }
  
  
  Z<-list(matrix(
    rep(c(delc,y_tilde,E,E_ind,Forc,M_at,M_up,M_lo,T_at,T_lo,Cum_D,Cum_E,
          Cum_dc,H,
          if(indic_stocks==1){c(Cum_div,pd,r.s,Cum_r.s,r.1)}
          ),
        nb.traj),
    model_sol$n.Z,nb.traj))
  
  Z[[nb.simul+1]]<-D                                                            #CHANGE if gamma0 not only on Tat
  Z[[nb.simul+2]]<-N
  if(indic_stocks==1){
    Z[[nb.simul+3]]<-D.div 
  }
  
  X<-list(model_sol$X)
  
  for (i in 2:nb.simul) {
    #Shocks
    Z[[nb.simul+1]][i,]    <-rgamma(nb.traj,rpois(nb.traj,param$ell0.D+
                                                    param$ell1.D[9]*Z[[i-1]][9,]),
                                    scale=param$mu_d)
    Z[[nb.simul+2]][i,]    <-rgamma(nb.traj,rpois(nb.traj,param$ell0.N+
                                                    param$ell1.N[9]*Z[[i-1]][9,]
                                                  +param$rho.N/param$mu_n*Z[[nb.simul+2]][i-1,])
                                    ,scale=param$mu_n)
    if(indic_stocks==1){
      Z[[nb.simul+3]][i,]    <-rgamma(nb.traj,rpois(nb.traj,param$ell0.div+
                                                      param$ell1.div[9]*Z[[i-1]][9,]),
                                      scale=param$mu_div) 
    }
    
    eta[1:model_sol$n.eta,]<-apply(eta,2,function(i)param$Phi%*%i)+
      rnorm(model_sol$n.eta*nb.traj)
    W                      <-rbind(eta,
                                   Z[[nb.simul+1]][i,],
                                   Z[[nb.simul+2]][i,],
                                   if(indic_stocks==1){
                                     Z[[nb.simul+3]][i,] 
                                   }
    )
    
    #Climate variables
    Z_1.0 <-apply(Z[[i-1]],2,
                  function(x){A1[[i-1]]%*%x+omega0[[i-1]]})
    Wz    <-apply(W,2,function(x){omega[[i-1]]%*%x})
    Z[[i]]<-Z_1.0+Wz
    X[[i]]<-rbind(Z[[i]],W)
  }
  delc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[1,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  y_tilde<-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[2,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  E      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[3,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  E_ind  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[4,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Forc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[5,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[6,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_up   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[7,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[8,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[9,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[10,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_D  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[11,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_E  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[12,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_dc <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[13,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  H      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[14,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  if(indic_stocks==1){
    Cum_div <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[15,]),1:nb.traj)),
                    ncol=nb.traj)[-1,]
    pd      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[16,]),1:nb.traj)),
                     ncol=nb.traj)[-1,]
    r.s     <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[17,]),1:nb.traj)),
                    ncol=nb.traj)[-1,]
    Cum_r.s <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[18,]),1:nb.traj)),
                    ncol=nb.traj)[-1,]
    r.1     <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[19,]),1:nb.traj)),
                    ncol=nb.traj)[-1,]
  }
  D      <-Z[[nb.simul+1]][-1,]
  N      <-Z[[nb.simul+2]][-1,]
  if(indic_stocks==1){
    D.div  <-Z[[nb.simul+3]][-1,] 
  }
  X      <-X[-1]
  
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=nb.simul-1)
  
  
  mylist<-c(list("delc"=delc,"y_tilde"=y_tilde,"E"=E,"E_ind"=E_ind,"Forc"=Forc,
               "M_at"=M_at,"M_up"=M_up,"M_lo"=M_lo,"T_at"=T_at,"T_lo"=T_lo,
               "Cum_D"=Cum_D,"Cum_E"=Cum_E,"Cum_dc"=Cum_dc,"H"=H),
               if(indic_stocks==1){list("Cum_div"=Cum_div,"pd"=pd,"r.s"=r.s,
                                        "Cum_r.s"=Cum_r.s,"r.1"=r.1)},
               list("D"=D,"N"=N),
               if(indic_stocks==1){list("D.div"=D.div)},
               list("X"=X))
  return(mylist)
}

################################################################################*LAPLACE TRANSFORMS FCTS
                                                                               #CHANGE if add gamma0
##Proposition 1: Laplace transform of W - Functions for 'a.w' and 'b.w'
#*U is a dim(W)*1 vector

a1.w.fct.1<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta],n.eta,1)
  u.d  <-U[n.eta+1]
  u.n  <-U[n.eta+2]
  if(indic_stocks==1){
    u.div<-U[n.eta+3] 
  }
  
  a.w<- t(u.eta)%*%u.eta/2+
    param$ell0.D*u.d*param$mu_d/(1-u.d*param$mu_d)+
    param$ell0.N*u.n*param$mu_n/(1-u.n*param$mu_n)+
    if(indic_stocks==1){
      param$ell0.div*u.div*param$mu_div/(1-u.div*param$mu_div)
    }else{
      0
    }
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    a.w[u.d>=1/param$mu_d]     <- NaN
    a.w[u.n>=1/param$mu_n]     <- NaN
    if(indic_stocks==1){
      a.w[u.div>=1/param$mu_div] <- NaN 
    }
  }
  
  
  return(a.w)
}

b1.w.fct.1<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta],n.eta,1)
  u.d  <-U[n.eta+1]
  u.n  <-U[n.eta+2]
  if(indic_stocks==1){
    u.div<-U[n.eta+3] 
  }
  
  a                                 <- matrix(0,model_sol$n.Z+model_sol$n.W,1)
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta)]        <- t(param$Phi)%*%u.eta
  b                                 <- param$ell1.D*
    u.d*param$mu_d/(1-u.d*param$mu_d)+
    model_sol$ell1.N.tilde*
    u.n*param$mu_n/(1-u.n*param$mu_n)+
    if(indic_stocks==1){
      param$ell1.div*
        u.div*param$mu_div/(1-u.div*param$mu_div) 
    }else{
      0
    }
  
  b.w <- a + b
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    b.w[u.d>=1/param$mu_d]     <- 100000
    b.w[u.n>=1/param$mu_n]     <- 100000
    if(indic_stocks==1){
      b.w[u.div>=1/param$mu_div] <- 100000 
    }
  }
  
  
  return(b.w)
}

################################################################################
                                                                               #CHANGE if add gamma0
#Proposition 1.N: Laplace transform of W - Functions for 'a.w' and 'b.w'
#*U is a dim(W)*N vector

a1.w.fct.N<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-U[n.eta+1,]
  u.n  <-U[n.eta+2,]
  if(indic_stocks==1){
    u.div<-U[n.eta+3,] 
  }
  
  a.w<-0
  a.w[1:dim(U)[2]]<- 0.5*matrix(apply(u.eta*u.eta,2,sum),
                                dim(U)[2],1)+
    matrix(param$ell0.D*u.d*param$mu_d/(1-u.d*param$mu_d),dim(U)[2],1)+
    matrix(param$ell0.N*u.n*param$mu_n/(1-u.n*param$mu_n),dim(U)[2],1)+
    if(indic_stocks==1){
      matrix(param$ell0.div*u.div*param$mu_div/(1-u.div*param$mu_div),dim(U)[2],1) 
    }else{
      0
    }
  a.w<-matrix(a.w,dim(U)[2],1)
  
  return(a.w)
}

b1.w.fct.N<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-matrix(U[n.eta+1,],ncol=dim(U)[2])
  u.n  <-matrix(U[n.eta+2,],ncol=dim(U)[2])
  if(indic_stocks==1){
    u.div<-matrix(U[n.eta+3,],ncol=dim(U)[2]) 
  }
  
  
  a                                      <- matrix(0,
                                                   model_sol$n.Z+model_sol$n.W,
                                                   dim(U)[2])
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta),]            <- t(param$Phi)%*%u.eta
  
  b                                      <- matrix(param$ell1.D%o%
                                                     (u.d*param$mu_d/
                                                        (1-u.d*param$mu_d))+
                                                     model_sol$ell1.N.tilde%o%
                                                     (u.n*param$mu_n/
                                                        (1-u.n*param$mu_n))+
                                                     if(indic_stocks==1){
                                                       param$ell1.div%o%
                                                         (u.div*param$mu_div/
                                                            (1-u.div*param$mu_div)) 
                                                     }else{
                                                       0
                                                     },
                                                   model_sol$n.Z+model_sol$n.W,
                                                   dim(U)[2])
  
  b.w <- a + b
  
  return(b.w)
}

################################################################################
##Proposition 2: One-Period ahead Laplace transform of X
#*A1, omega and omega0 start at t=1 (2020)
#*t=0 is for t=2015, t=1 is for t=2020, etc.
#*MAX=98, with 99=inf period
#*U is a matrix of dimension "model_sol$n.X*1" 

a1.fct.1<-function(model_sol,U,t){
  omega.0   <-model_sol$omega0
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z],model_sol$n.Z,1)
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)],model_sol$n.W,1)
  
  a <-t(Uz)%*%omega.0[[t+1]]+a1.w.fct.1(model_sol,Uw+t(omega[[t+1]])%*%Uz)
  
  return(a)
}

b1.fct.1<-function(model_sol,U,t){
  A.1       <-model_sol$A1
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z],model_sol$n.Z,1)
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)],model_sol$n.W,1)
  
  b <-rbind(t(A.1[[t+1]])%*%Uz,matrix(0,model_sol$n.W,1))+
      b1.w.fct.1(model_sol,Uw+t(omega[[t+1]])%*%Uz)
  
  return(b)
}

################################################################################
##Proposition 2.N: One-Period ahead Laplace transform of X
#U is a model_sol$n.X*N matrix

a1.fct.N<-function(model_sol,U,t){
  omega.0   <-model_sol$omega0
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],
             model_sol$n.W,dim(U)[2])
  
  a <-t(Uz)%*%omega.0[[t+1]]+a1.w.fct.N(model_sol,Uw+t(omega[[t+1]])%*%Uz)
  
  
  return(a)
}

b1.fct.N<-function(model_sol,U,t){
  A.1       <-model_sol$A1
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):model_sol$n.X,],model_sol$n.W,dim(U)[2])
  
  b <-rbind(t(A.1[[t+1]])%*%Uz,matrix(0,model_sol$n.W,dim(U)[2]))+
      b1.w.fct.N(model_sol,Uw+t(omega[[t+1]])%*%Uz)
  
  return(b)
}

################################################################################*INFINITE AND GN
##Functions for infinite a,b and Gauss-Newton
#*U can be of dimension dim(U)=X*N

a1.fct.inf<-function(model_sol,U){
  omega.0.inf  <-model_sol$omega0.inf
  omega.inf    <-model_sol$omega.inf
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],model_sol$n.W,
             dim(U)[2])
  
  a.inf<-t(Uz)%*%omega.0.inf+a1.w.fct.N(model_sol,Uw+t(omega.inf)%*%Uz)
  
  return(a.inf)
}

b1.fct.inf<-function(model_sol,U){
  omega.inf  <-model_sol$omega.inf
  A.1.inf    <-model_sol$A1.inf
  
  Uz<-matrix(U[1:model_sol$n.Z],model_sol$n.Z,1)
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)],model_sol$n.W,1)
  
  b.inf<-rbind(t(A.1.inf)%*%Uz,matrix(0,model_sol$n.W,1))+
    b1.w.fct.1(model_sol,Uw+t(omega.inf)%*%Uz)
  
  return(b.inf)
}

#######################################
##Gauss Newton mu_u1
#*x0 is a matrix of dimension "model_sol$n.X*1", ini guess
#*return nb of iterations, dev, mu_u1.inf, mu_c1, list dev

GN.function<-function(model_sol,x0){
  param     <-model_sol$parameters
  A1.inf    <-model_sol$A1.inf
  omega.inf <-model_sol$omega.inf
  #Jacobian
  e1.X    <-matrix(0,model_sol$n.Z+model_sol$n.W,1)
  e1.X[1] <-1
  mu_c1   <-e1.X
  mu_u1   <-x0
  eps     <-param$eps.GN
  J       <-matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z+model_sol$n.W)
  dev     <-matrix(1,model_sol$n.Z+model_sol$n.W,1)
  ite     <-0
  listdev <-list(dev)
  tol     <-max(abs(dev*param$tol.GN))
  while((max(abs(dev))>tol)&(ite<50)){
    b  <- b1.fct.inf(model_sol,(1-param$gamma)*(mu_u1+mu_c1))
    dev<--mu_u1+param$delta/(1-param$gamma)*b
    f  <-dev
    for (j in 1:(model_sol$n.Z+model_sol$n.W)){
      mu_u1.perturb   <- mu_u1
      mu_u1.perturb[j]<- mu_u1[j]+eps
      U.perturb       <- (1-param$gamma)*(mu_u1.perturb+mu_c1)
      b.perturb       <- b1.fct.inf(model_sol,U.perturb)
      dev.perturb     <--mu_u1.perturb+param$delta/(1-param$gamma)*b.perturb
      f.prime         <- dev.perturb
      J[,j]           <-(f.prime-f)/eps
    }
    mu_u1         <-mu_u1-ginv(J)%*%dev                                         #generalized inverse
    ite           <-ite+1
    listdev[[ite]]<-dev
  }
  mylist<-list("ite"=ite,"dev"=dev,"mu_u1"=mu_u1,"listdev"=listdev,
               "mu_c1"=mu_c1)
  return(mylist)
}

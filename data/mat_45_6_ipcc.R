#Scenarios to compute m_0 (RF equation)

data6<-read.table("../Data/State_Vector_Data/RCP6_MIDYR_CONC_R.txt")

mat6<-data6[-(1:3),1:2]

colnames(mat6)<-c("year","CO2eq_ppm")

#ppm into GtC:

mat6[,2]<-as.numeric(mat6[,2])*2.123

plot(mat6[,1],mat6[,2],type="l")
abline(v=2100)
#########################
data45<-read.table("../Data/State_Vector_Data/RCP45_MIDYR_CONC_R.txt")

mat45<-data45[-(1:3),1:2]

colnames(mat45)<-c("year","CO2eq_ppm")

#ppm into GtC:

mat45[,2]<-as.numeric(mat45[,2])*2.123

plot(mat45[,1],mat45[,2],type="l")
abline(v=2100)


submat45<-mat45[256:336,]
submat6<-mat6[256:336,]

plot(submat45[,1],submat45[,2],type="l",ylim=c(range(submat45[,2],submat6[,2])))
lines(submat6[,1],submat6[,2],col="red")

mean45<-mean(submat45[,2])
mean6<-mean(submat6[,2])

mean(c(mean45,mean6))
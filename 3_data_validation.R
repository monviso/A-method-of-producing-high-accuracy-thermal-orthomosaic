####VALIDATION

##check predicted values
Tk<-raster::extract(n,log_sf_val,buffer=0.1,fun=mean,na.rm=F)#Extract Tk with validation loggers, 0.1m buffer
comp2<-Tk
comp2<-data.frame(V1=comp2)

#Produce a table (comp_long) with AB, MAB, MAB sd, RMSE, r
comp_long<-gather(comp2,d,predT,V1)  
comp_long$logger<-rep(log14_sf$logger[c(valid)],1) 
comp_long$meanT<-rep(log14_sf$temp[c(valid)],1)
comp_long$diff<-comp_long$predT-comp_long$meanT
comp_long$MAB<-rep(aggregate(diff~d,comp_long,FUN=mean)[,2],each=1)
f1<-function(x) (sd = sd(x))
comp_long$MAEsd<-rep(aggregate(diff~d,comp_long,FUN=f1)[,2],each=1)


library(ModelMetrics)
library(dplyr)
g<-comp_long[-9,]%>%
  group_by(d) %>%
  summarise(
    RMSE = rmse(meanT, predT)
    ,R = cor(meanT, predT)
  )
comp_long $RMSE<-rep(g$RMSE,each=1)
comp_long $R<-rep(g$R,each=1)


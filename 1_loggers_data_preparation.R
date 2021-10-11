##load libraries
library(sf)
library(raster)
library(ggplot2)
library(tidyr)

##preparing sf object with loggers position
log_14<-read.table("~/jan20a.txt",sep=",",
                   colClasses = c(rep("NULL",6),"character",rep("NULL",3),"character","numeric","numeric"),header=T)
log14_sf<-st_as_sf(log_14,coords = c("long","lat"))
log14_sf<-log14_sf[log14_sf$logger!="T" & log14_sf$logger!="TK" & log14_sf$logger!="GCP",] ##removing from table unwanted point (Ground control points and take-off)
st_crs(log14_sf)<-4326 #set projection WGS84

##preparing water Tk recored by loggers during the flight
hotsp_long<-read.csv(file="~/hotsp_long.csv",sep=",",header = T)
hotsp_long$date.time<-as.POSIXct(hotsp_long$date.time,format="%Y-%m-%d %H:%M:%S")
temp_14<-subset.data.frame(hotsp_long,date.time>"2020-01-22 07:03:00" & date.time<"2020-01-22 07:10:00")#flight time
temp_14$temp<-as.numeric(temp_14$temp)

##check Tk trend during the flight 
ggplot(data= subset.data.frame(temp_14,pos!="air1"&pos!="air2"),aes(x=date.time,y=temp,color=pos))+
  geom_line()+
  labs(x="hour (AM)",y="temperature (°C)",color="loggers ID",title = "Survey 1")+
  theme_bw()
mm<-function(x){max(x)-min(x)}
d<-aggregate(temp~pos,temp_14,FUN=mm)#table of difference between max and min Tk recorded during the flight
##IMPORTANT: check for eventual water Tk variation during flight time before proceeding with the use of mean Tk  

##adding mean temperature during the flight time of each logger
mn<-aggregate(temp~pos,temp_14,FUN=mean)
log14_sf<-merge(log14_sf,mn,by.x="logger",by.y="pos")

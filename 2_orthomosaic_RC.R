
#GRAYSCALE ORTHOMOSAIC ELABORATION
###importing original grayscale orthomosaic

p<-("~/orthomosaic/")
pict_name<-"gray_ICE_rect1.tif"#use the 15-35 range orthomosaic
file<-paste0(p,"gray/",pict_name)
img<-raster(file)
img<-img/255#0-1 range
water<-read_sf("~/water_mask/water_mask_jan20.shp")#importing water mask; to be manually created in GIS
img2<-mask(img, water)
plot(img2)

#extract for each logger the original grayscale value (GVf, equation 4)
GVf<-raster::extract(img2,log14_sf,buffer=0.25,fun=mean,na.rm=T)#buffer of 0.25 m
log14_sf$GVf<-GVf

#importing RC equations built in laboratory
rgbmodel<-read.csv("~/rgbmodel_gray.csv")
rgbmodel$range<-as.factor(rgbmodel$range)
levels(rgbmodel$range)<-c("15-35 °C","15-45 °C","20-45 °C","25-45 °C" )

#build RC1535 
lin_1535<-lm(bath~gray,data = subset.data.frame(rgbmodel,range=="15-35 °C"))

#build the revers of RC1535
lin_1535rev<-lm(gray~bath,data = subset.data.frame(rgbmodel,range=="15-35 °C"))

#compute the expected grayscale value as for laboratory condition for the mean Tk recored by loggers (GVlab, equation 4)
val<-data.frame(bath=log14_sf$temp)
log14_sf$GVlab<-predict(lin_1535rev,val)

#Compute Grayscale Correction Values
log14_sf$GCVs<-log14_sf$GVlab-log14_sf$GVf

##select computation and validation loggers
valid<-c(2,6,8,11,13,14,17,21,23,26) #loggers selected to be homogenously distributed in the reach
log_sf_val<-log14_sf[valid,] #validation loggers
log_sf_com<-log14_sf[-valid,] #computation loggers

####QUESTION: 
#GCVs can be considered homogeneous?#####
##to consider: for the range 15-35°C(amplitude 20°C) a variation of 0.1 of grayscale intesity 
#correspond to a variation of ~2.6 °C as for RC1535

ggplot(log14_sf,aes(x=logger,y=GCVs))+
  geom_point(group=1)
sd(GCVs)

###ANSWER

#YES: use a single GCV to correct the Gf, i.e. mean GCV
GCV<-mean(log_sf_com$GCVs)

#NO: determine GCVs spatial distribution using GAM
library(mgcv)
##hot to use a GAM with soap film smoother refr to:
#set boundaries of the water area
g<-st_coordinates(water)
g<-as.data.frame(g)
g_sf<-st_as_sf(g,coords=c("X","Y"))
st_crs(g_sf)<-4326

#the GCV values to be assumed by the boundary is the GCV of the closest logger
min_dist<-st_distance(g_sf,log14_sf)

ddd<-apply(min_dist, 1, FUN=min)
r<-c()
for (i in 1:nrow(min_dist)) {
  s<-as.numeric(min_dist[i,])
  d<-which(s==ddd[i])
  r[i]<-as.numeric(d)
  
}
bd_gcv<-c()
for (i in 1:length(r)) {
  w<-log14_sf$diff[which(log14_sf$seq==r[i])]
  bd_gcv[i]<-w
}

##setting boundary for the smoother
bound <- list(list(long = g[,1], lat = g[,2], f=bd_gcv))##f=the GCV value of the closest logger
N <- 9 #number of knots
gx <- seq(min(g[,1]), max(g[,1]), len = N)
gy <- seq(min(g[,2]), max(g[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("long","lat")
gp_sf<-st_as_sf(gp,coords = c("long","lat"))
st_crs(gp_sf)<-4326
knots <- gp[with(gp, inSide(bound, long, lat)), ]
out<-autocruncher(bound,knots,x="long",y="lat") #download autocruncher form:
knots_sf<-st_as_sf(knots,coords = c("long","lat"))
st_crs(knots_sf)<-4326

#plot of knots, only the knots inside the boundary will be used
ggplot()+
  geom_sf(data=water)+
  geom_sf(data=gp_sf,color="blue")



#latitude and longitude of the computation loggers to be used as covariates
crd<-st_coordinates(log_sf_com)
log_sf_com$long<-crd[,1]
log_sf_com$lat<-crd[,2]

##computing equation5
gam1<-gam(GCVs~s(long,lat,bs="so",xt=list(bnd=bound))+s(GVf,k=3),data=log_sf_com ,knots = knots[-c(out),] ,method = "REML")
summary(gam1)
gam.check(gam1)

#predict GCVs distribution and values using equation5 (gam1)
pre_img<-data.frame(long=img_df$x,lat=img_df$y,GVf=img_df$gray_ICE_rect1)#covariates form the raster
pred<-predict.gam(gam1, newdata =pre_img,type="response")#prediction

#shape back the predicted GCVs to a raster of GCVs
pred_rs<-raster(matrix(pred,ncol=ncol(img2),byrow = T))
extent(pred_rs)<-extent(img2)
crs(pred_rs)<-4326
plot(pred_rs)


##PERFORM RASTER CORRECTION:

#Mean GCV
##compute raster with mean GCV
bw.values<-data.frame(gray=img_df$gray_ICE_rect1+ GCV) ##add mean GCV to each pixel's GVf
pred2<-predict(lin_1535,bw.values)#predict Tk values

#shape back the predicted Tk to a raster 
pred<-matrix(pred2,nrow=nrow(img2),byrow = T)
n<-raster(pred)
crs(n)<-4326
extent(n)<-extent(img2)
Tk_raster_meanGCV<-mask(n,water)#final Tk raster using a single GCV


#GAM modeled GCVs
img_sum<-img2+pred_rs#adding original raster (GVf) with GCVs raster

#predict Tk values
img_sum_df<-as.data.frame(img_sum,xy=T)
bw.values<-data.frame(gray=img_sum_df$layer)
pred2<-predict(lin_1535,bw.values)

#shape back the predicted Tk to a raster 
pred<-matrix(pred2,nrow=nrow(img2),byrow = T)
n<-raster(pred)
crs(n)<-4326
extent(n)<-extent(img2)
Tk_raster_GAMGCVs<-mask(n,water)

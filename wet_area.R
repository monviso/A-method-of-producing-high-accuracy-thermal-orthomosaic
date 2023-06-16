wet_area<-function(wt_path,rgb_path,thermal_path,out_path,crop_path,NAval=0,classID=c("w","t"),ref_crs=4326,lat_long=T){
  if(require(terra)==F|require(sf)==F|require(mgcv)==F|require(stars)==F){stop("install required packages (terra, sf, mgcv, stars)!")}else{
p<-read_sf(wt_path)
p$wt01<-ifelse(p$wt==classID[1],1,0) #assign 0-1 code
rgb<-rast(rgb_path)
thermal<-rast(thermal_path)
cr<-read_sf(crop_path)
NAflag(rgb)<-NAval
NAflag(thermal)<-NAval
rgb<-crop(rgb,cr)
thermal<-crop(thermal,cr)

rgb<-resample(rgb,thermal)
ml<-c(rgb,thermal)
ml<-mask(ml,cr)

s12_df<-data.frame()
for (i in 1:nrow(p)) {
  #show(i)
  ms<-p[i,]
  rs<-crop(ml,ms)
  rs<-mask(rs,ms)
  s_df<-as.data.frame(rs,xy=T,na.rm=F)
  na<-which(is.na(s_df[,3])&is.na(s_df[,4])&is.na(s_df[,5])&is.na(s_df[,6]))
  if(length(na)==0|length(na)==nrow(s_df)){s_df<-s_df}else{s_df<-s_df[-na,]}
  
  #z<-which(s_df[,3]==0|s_df[,6]==0)
  #if(length(z)==0|length(z)==nrow(s_df)){s_df<-s_df}else{s_df<-s_df[-z,]}
  s_df$wt01<-ifelse(ms$wt01==1,1,0)
  s12_df<-rbind(s12_df,s_df)
}

s12_df<-st_as_sf(s12_df,coords = c("x","y"),remove=F)
s12_df$row<-seq(1,nrow(s12_df))
colnames(s12_df)<-c("long","lat","R","G","B","th","wt01","geometry","row")
s12_df$wt01<-as.factor(s12_df$wt01)
type<-c(levels(s12_df$wt01))
p_tr<-st_sf(st_sfc())
p_vl<-st_sf(st_sfc())
st_crs(p_tr)<-ref_crs
st_crs(p_vl)<-ref_crs

for (i in 1:length(type)) {
  df<-s12_df[which(s12_df$wt01==type[i]),]
  df_tr<-df[sample(nrow(df), nrow(df)*(2/3)), ]
  df_vl<-df[-which(df$row %in% df_tr$row),]
  p_tr<-rbind(p_tr,df_tr)
  p_vl<-rbind(p_vl,df_vl)
}

if(lat_long==T){
gam_wt_CAT<-bam(wt01~s(th)+s(R,G,B)+te(long,lat),data=p_tr,family="binomial",nthreads = 4,discrete = T)
}else{
gam_wt_CAT<-bam(wt01~s(th)+s(R,G,B),data=p_tr,family="binomial",nthreads = 4,discrete = T)
  
}

r01df<-as.data.frame(ml,xy=T,na.rm=F)#transfm it into a daframe format 
r01df$row<-rownames(r01df)
colnames(r01df)<-c("long","lat","R","G","B","th","row")

r01df$wt_prd<-predict.gam(gam_wt_CAT,r01df,type = "response")#predict class of each pixel
r01df$wt_prd01<-ifelse(r01df$wt_prd>=0 & r01df$wt_prd<0.5,0,1)
dwa_m<-matrix(r01df$wt_prd01,ncol=ncol(ml),byrow=T)
dwa<-rast(dwa_m)
ext(dwa)<-ext(ml)
crs(dwa)<-crs(ml)

library(stars)
wetted <- st_as_stars(dwa)%>%st_as_sf(merge = TRUE)#produce an shp
wetted$area<-st_area(wetted)

maina<-wetted[wetted$lyr.1==1,]
main<-maina[which.max(maina$area),]

write_sf(main, paste0(out_path,"wetted_polygon.shp"))
writeRaster(dwa,paste0(out_path,"wetted_polygon.shp"),overwrite=T)

pw<-p_vl[which(p_vl$wt01==1),]
pt<-p_vl[which(p_vl$wt01==0),]
w<-terra::extract(dwa,pw)
t<-terra::extract(dwa,pt)
#compute True Positivives, False negative, False Positives

Tp<-nrow(w[which(w$lyr.1==1),])
Fn<-nrow(w[which(w$lyr.1==0),])+nrow(w[is.na(w$lyr.1),])
Fp<-nrow(t[which(t$lyr.1==1),])
#Tp;Fn;Fp
P<-Tp/(Tp+Fp)#precison
R<-Tp/(Tp+Fn)#recall

F1<-2*((P*R)/(P+R))#F1
#P;R;F1

perf<-data.frame(tr_pixel=nrow(p_tr), val_pixel=nrow(p_vl), Tp=Tp, Fn=Fn,Fp=Fp,P=P,R=R, F1=F1)
write.csv(perf,paste0(out_path,"class_performance.csv"))

  }
}
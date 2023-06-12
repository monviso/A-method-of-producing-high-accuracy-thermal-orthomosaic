UAVtoraw<-function(pt,rgb=T,gamma=F,pt.vign.corr,pt.dest,pt.exif){
  if(require(terra)==F){stop("install package terra first!")}else{
    
    g<-list.files(pt)
    vign_corr<-terra::rast(pt.vign.corr)
    for (i in 1:length(g)){
      suppressWarnings({
        show(g[i])
        r<-terra::rast(paste0(pt,g[i]))
        if(rgb=T & gamma=F){
          r2<-as.array(r/256)
          bw.values <- r2[,,1]*0.2162 + r2[,,2]*0.7152 + r2[,,3]*0.0722
          bw.values<-terra::rast(bw.values)
          terra::ext(bw.values)<-terra::ext(r)
          bw.values<-bw.values+vign_corr
        } else if(rgb=T & gamma=T){
          r2<-as.array(r/256)
          r2<-ifelse(r2<0.04045,r2/12.92,(r2+0.055)/1.055)
          bw.values <- r2[,,1]*0.2162 + r2[,,2]*0.7152 + r2[,,3]*0.0722
          bw.values<-ifelse(bw.values<=0.0031308,bw.values*12.92,(1.055*bw.values^(1/2.4))-0.055)
          bw.values<-terra::rast(bw.values)
          terra::ext(bw.values)<-terra::ext(r)
          bw.values<-bw.values+vign_corr
          
        } else{bw.values<-r+vign_corr}
        
        terra::writeRaster(bw.values,filename =paste0(pt.dest,g[i],".tiff"),overwrite=T)
        
        meta<-system2(pt.exif,paste("-s","-g1",paste0(pt,g[i])),stderr = T)
        meta<-str_replace_all(meta,fixed(" "), "")
        meta<-gsub('["]','',meta)
        meta.df<-data.frame(meta)
        g2<-meta.df %>% tidyr::separate(meta,c("A","B"),sep = ":",extra = "merge")
        g2<-g2[-c(which(is.na(g2$B))),]
        row.names(g2)<-seq(1:length(g2$A))
        nec<-g2[c(20,21,28,33,34,51,52,60:69,89:91),]#number customized for DJI MDE (Lepton3.0)
        att<-c(paste0("-",nec$A,"=",nec$B),"-overwrite_original")
        system2(pt.exif,arg= paste0(att," ",paste0(pt.dest,g[i])))
        
      }) }
    
  }
  
}


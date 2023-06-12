vign_matrix<-function(pt,dest_pt=pt){
  if(require(terra)==F){stop("install package terra first!")}else{
    
  g<-list.files(pt)
  for (i in 1:length(g)){
    suppressWarnings({
    m<-terra::rast(paste0(pt,g[i]))
    if(dim(m)[3]>1){stop("images must be single layer (e.g.: raw, grayscale)")}
    else if(i==1){m2<-m}else{m2<-(m2+m)/i}
    }) }
  mm<-matrix(m2,ncol = ncol(m),byrow = T)
  
  mx<-max(mm,na.rm=T)
  mc<-mx-mm
  cr<-terra::rast(mc)
  terra::ext(cr)<-terra::ext(m)
   
  terra::writeRaster(cr,paste0(dest_pt,"vign_matrix.tiff"))
    
  
}
}


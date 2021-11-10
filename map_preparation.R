## necessary modules: module load gcc/6.3.0; module load GEOS/3.8.1 UDUNITS/2.2.26 R/3.5.0; module load PROJ/7.0.1; module load GDAL/2.4.2; LIBRARY_PATH=/apps/GCC/6.3.0/lib64:$LIBRARY_PATH
# R --vanilla
library(gstat)
library(sp)
library(raster)
library(spData)
library(dplyr)
library(sf)
library(tmap)
library(tmaptools)
library(maptools)
library(automap)
library(spatialEco)
library(splancs)
library(ggplot2)
library(rgdal)
library(mapview)
options("scipen"=100)
Sys.setenv("PROJ_LIB" = "/apps/PROJ/7.0.1/share/proj/")
'%ni%' <- Negate('%in%')

### spatial map preparation
load(file="data/colnams6.Robject")
load(file=paste("data/coln",sep=""))
koord<-list()
for (j in (1:length(coln))) { koord[[j]]<-snam4[which(snam4[,2]==coln[j]),9:10][1,]  }
koord<-matrix(as.numeric(unlist(do.call(rbind,koord))),ncol=2)
koord2<-paste(koord[,1],", ",koord[,2],", ",0,sep="")

# map of Africa
ext<-c(-17,35,-12,17)
africa = world %>% filter(continent == "Africa" & name_long%in% c("Tanzania","Democratic Republic of the Congo","Kenya","Chad","Senegal","Mali","Mauritania","Benin","Niger","Nigeria","Cameroon","Togo","Ghana","Côte d'Ivoire","Guinea","Guinea-Bissau","Liberia","Sierra Leone","Burkina Faso","Central African Republic","Republic of the Congo","Gabon","Equatorial Guinea","Zambia","Angola","Burundi","The Gambia","Uganda","Rwanda","South Sudan","Sudan","Ethiopia"))
africa<-crop_shape(africa, extent(ext), polygon = FALSE)
crs    <- CRS("+proj=longlat +ellps=WGS84") 

## coordinates
koord3<-SpatialPoints(coords= koord, proj4string = crs)

## chimps area (map from IUCN)
area.c <- readOGR("Pan_trogl_range_2016.shp")
area.d.points <- fortify(area.c[1,])

## extend the area
mm=0
amp<-list()
for (j in 1:length(area.c@polygons)) {
  mp<-area.c@polygons[[j]]@Polygons
  for (i in (1:length(mp))) { 
    mm=mm+1
    np<-list(Polygon(as.data.frame(mp[[i]]@coords)))
    np<-SpatialPolygons(list(Polygons(np,ID=1)))
    np@bbox<-t(matrix(ext,2,2));np = st_buffer(st_as_sf(np), 1.5)
    np<-as_Spatial(np$geometry,ID=paste("no_",mm,sep=""))
    np@bbox<-t(matrix(ext,2,2))
    amp[[mm]]<-np
  }
}

np2<-unionSpatialPolygons(SpatialPolygons(lapply(amp, function(x){x@polygons[[1]]})),IDs=paste("no_",c(1:length(amp)),sep=""))
alco<-list()
for (i in (1:length(amp))) { alco[[i]]<-SpatialPoints(coords= np2@polygons[[i]]@Polygons[[1]]@coords, proj4string = crs);  alco[[i]]@bbox<-t(matrix(ext,2,2)) }

## extended area to predict in
mp<-list(Polygon(as.data.frame(koord3)));names(mp)<-"p1"
ma.poly<-SpatialPolygons(list(Polygons(mp,ID="p1")))
ma.poly@bbox<-t(matrix(ext,2,2));ma.poly = st_buffer(st_as_sf(ma.poly), 4)
ma.poly<-as_Spatial(ma.poly$geometry)
ma.poly@bbox<-t(matrix(ext,2,2))
koord4<-SpatialPoints(coords= ma.poly@polygons[[1]]@Polygons[[1]]@coords, proj4string = crs)
rownames(koord4@bbox)<-rownames(koord3@bbox)

# extended grid
P <- SpatialPointsDataFrame(coords= koord3, data= as.data.frame(rep(0,length(koord3))),  proj4string = crs,bb=koord4@bbox)
P@bbox<-koord4@bbox
fgrd              <- as.data.frame(spsample(P, "regular", n=100000))

## loop for inout
ino<-list()
for (i in (1:length(alco))) { ino[[i]]<-inout(fgrd,as.data.frame(alco[[i]]),bound=T)}
ino<-do.call(cbind,ino)
ino<-ifelse(rowSums(ino)>0,T,F)
grd<-fgrd[ino,]
names(grd)       <- c("X", "Y");coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE
proj4string(grd) <- crs

save(grd,crs,fgrd,P,koord,koord3,koord2,koord4,africa,ext,coln,snam4,file=paste("data/spatial_stuff",sep=""))


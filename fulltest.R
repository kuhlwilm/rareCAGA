#!/usr/bin/r
print(Sys.time())
ip=unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" "))
print(ip)

## input: vcf file (bcftools suitable, "chr21" chromosome name convention)
infile=ip[1];labl=ip[2]

# if necessary, create results directory
dir.create("results")

# load the general functions and options
bps<-c("A","C","G","T")
options("scipen"=100)
'%ni%' <- Negate('%in%')
fifu<-function(input) { if(input!=".") { return(min(as.numeric(unlist(strsplit(input,split=","))))) } else { NA } }
yle=seq(0,49000000,1000000)
load(file="data/spatial_stuff")
yvec<-c(15:49)

# load the segregating positions
load(file="data/aprivatT9")

# sample names
snams<-system(paste("bcftools view -h ",infile," | tail -n 1",sep=""),intern=T)
snams<-unlist(strsplit(snams,split="\t"))[-c(1:9)]


############################################################################################################
## step one, filtering and overlap with private alleles
############################################################################################################

alovr<-list()
save(alovr,file=paste("results/",labl,sep=""))


## this is iterating in 1Mbp steps to get the overlap with ~private SNPs
## stepwise procedure to avoid memory failure
for (yy in yvec) {
  print(yy)
  # load genotypes and perform basic filtering
    agt<-system(paste("bcftools view -Ou -M2 -R data/derived_9sites",yy,".bed ",infile," | bcftools query -u -f '%POS %REF %ALT [%GT,] [%DP,] [%GQ,] [%AD-]\n' -",sep=""),intern=T)
    print("loaded")
    if (length(agt)==0) {  next } 
      agt<-do.call(rbind,strsplit(agt,split=" "))
      nfo<-agt[,1:3,drop=F]
      gt<-do.call(rbind,strsplit(agt[,4],split=","))
      dp<-matrix(as.numeric(do.call(rbind,strsplit(agt[,5],split=","))),nrow=nrow(nfo))
      gq<-matrix(as.numeric(do.call(rbind,strsplit(agt[,6],split=","))),nrow=nrow(nfo))
      ad<-matrix(do.call(rbind,strsplit(agt[,7],split="-")),nrow=nrow(nfo))
      ad<-matrix(as.numeric(unlist(apply(ad,c(1:2),fifu))),nrow=nrow(ad))
      
      # filtering for coverage 2-99 and genotype quality >=20
      GQ<-gt
      GQ[which(dp<1)]<-"./."
      GQ[which(dp>99)]<-"./."
      GQ[which(gq<20)]<-"./."
      # heterozygous sites can be filtered for allele balance if such information is available
      GQ[which(gt%in%c("0/1","1/2","2/1","0/2","0/3","1/3","3/1","3/2","2/3","0/4","1/4","2/4","3/4","4/1","4/2","4/3") & (((ad/dp)<0.1)|ad==1))]<-"./."
      gtQ<-GQ    
      
      # make sure only SNPs are used
      selc<-which(nfo[,2]%in%bps & nfo[,3]%in%bps)
      nfo2<-nfo[selc,,drop=F]
      gtAQ<-gtQ[selc,,drop=F]
      # transform genotypes to simple format
      gtAQ[which(gtAQ=="1/1")]<-2;gtAQ[which(gtAQ=="0/1")]<-1;gtAQ[which(gtAQ=="0/0")]<-0;gtAQ[which(gtAQ=="./.")]<-NA
      print("filtered")  
      
      ## merge with private dataset
      # if necessary, remove sites that are monomorphic
      asi<-asites[[yy]]
      ssi<-merge(asi[,-2],nfo2[,-2,drop=F],by.x=1,by.y=1,all=T)
      remo<-ssi[which(ssi[,3]!="." & unlist(as.character(ssi[,2]))!=ssi[,3]),1]
      privat<-apri[[yy]]
      if (length(remo)>0) { gtAQ<-gtAQ[-which(nfo2[,1]%in%remo),,drop=F];nfo2<-nfo2[-which(nfo2[,1]%in%remo),,drop=F];privat<-privat[-which(asites[[yy]][,1]%in%remo),,drop=F];asi<-asi[-which(asites[[yy]][,1]%in%remo),,drop=F] }
      selsi<-which(asi[,1]%in% nfo2[,1])
      if(length(selsi)==0) { next }
      # get the subset of informative sites in the data
      segs2<-gtAQ[match(asi[,1],nfo2[,1]),,drop=F]
      segs2a<-matrix(as.numeric(segs2),nrow=nrow(segs2))
      for (i in (1:ncol(segs2))) {     segs2a[,i]<-ifelse(segs2a[,i]%in%c(1,2),1,segs2a[,i]) }
      # calculate the statistics
      # 1) Observed sites in both; 2) Observed & derived in Reference; 3) Observed & derived in Test; 4) Identical Ref/Test; 5) Total derived in Ref; 6) Derived in both
      ovrla<-list()
      for (i in (1:ncol(segs2a))) {
        ovrla[[i]]<-list()
        for (j in (1:ncol(privat))) {
          sgs<-cbind(privat[,j],segs2a[,i])
          ovrla[[i]][[j]]<-c(length(which(!is.na(rowSums(sgs))==T)),length(which(!is.na(rowSums(sgs))==T& sgs[,1]==1)),length(which(!is.na(rowSums(sgs))==T& sgs[,2]==1)),length(which(sgs[,1]==sgs[,2])),length(which(sgs[,1]==1)),length(which(sgs[,1]==sgs[,2] & sgs[,1]==1)))
        }
        ovrla[[i]]<-do.call(rbind,ovrla[[i]])
      }
      alovr[[yy]]<-ovrla  
      save(alovr,file=paste("results/",labl,sep=""))
      print("done")
}

print(Sys.time())
#q()

############################################################################################################
## step two, aggregate data and apply spatial model
############################################################################################################

## if successful, generate the full table
print("full table")
load(file=paste("results/",labl,sep=""))
y=0; repeat { y=y+1;if (y==49) { print("no data");q() }; if(length(alovr[[y]])==0) { next } else {strt=y;break}  }
a1<-alovr[[strt]]
  for (i in (strt+1):(length(alovr))) {
  for (j in (1:length(a1))) {   if(length(alovr[[i]][[j]])>0) { a1[[j]]<-a1[[j]]+alovr[[i]][[j]] } }
}
ovrlai<-list()
for (j in (1:length(a1))) {
  ovrlai[[j]]<-a1[[j]][,6]/a1[[j]][,5]
}
ovrlai<-do.call(cbind,ovrlai)
tmpa<-ovrlai
tmpa[is.na(tmpa)]<-0
## this will save the aggregated table that can be used for downstream analysis
save(tmpa,ovrlai,a1,file=paste("results/R1_",labl,sep=""))

## spatial analysis
# needs a bunch of R libraries for kriging and creating the maps
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
library(TSCS)
#Sys.setenv("PROJ_LIB" = "/apps/PROJ/7.0.1/share/proj/")

# load the data
load(file="data/spatial_stuff")
load(paste("results/R1_",labl,sep=""))

# loop for kriging and map creation per sample
print("maps")
allTs<-list()
for (i in (1:ncol(tmpa))) {
  print(i)
  # take the data for one sample
  subs<-as.data.frame(as.vector(unlist(tmpa[,i])))
  if (sum(subs)==0) { next }
  # mark the best fitting sampling location and create a spatial object for all locations
  mc<-rep("grey40",nrow(tmpa));mc[which(subs==max(subs))]<-"black"
  P <- SpatialPointsDataFrame(coords= koord3, data= subs,  proj4string = crs)
  P$colr=mc;P@bbox<-t(matrix(ext,2,2))
  P$X=koord[,1];P$Y=koord[,2]
  
  # fit variogram and perform kriging on the grid
  f.1 <- as.formula(unlist(subs) ~ X + Y) 
  var.smpl <- variogram(f.1, P, data=P , cloud = FALSE)
  dat.fit  <- fit.variogram(var.smpl, fit.ranges = F, fit.sills = F, vgm( model="Sph"))
  dat.krg <- krige( f.1, P, grd, dat.fit )
  # which is the best spatial fit?
  cetr<-SpatialPoints(dat.krg@coords[which(dat.krg$var1.pred==max(dat.krg$var1.pred)),,drop=F],  proj4string = crs,bbox=koord4@bbox);cetr$colr="red"
  
  # rasterize this object
  r <- raster(dat.krg,values=T)
  # put it on the african map
  r.m <- mask(r, africa);  r.m<-extend(r.m,extent(ext))
  # make sure data is there
  r.m<-readAll(r.m)
  
  # the plot: krige raster, country borders, sampling locations, best fit
  t1<-tm_shape(r.m) +
    tm_raster(n=15, palette="RdBu",title="Similarity",midpoint=NA,style="cont") +
    tm_shape(africa) + 
    tm_borders() +
    tm_shape(P) + tm_dots(size=0.2,col="colr") +
    tm_legend(legend.outside=TRUE) +
    tm_layout(title=paste(snams[i],sep="")) +     tm_shape(cetr) + tm_dots(size=0.8,shape=3,col="red") 

  # collect plots for all samples
  allTs[[i]]<-t1
}


## plot 10 samples per page
pm<-ceiling(ncol(tmpa)/10)
allT3<-list();va=0
for (l in (1:pm)) {
  allT3[[l]]<-list();k=0
  for (j in 1:10) { va=va+1; if (va>ncol(tmpa)) { break }; if (length(allTs[[va]])==0) { next };      k=k+1;allT3[[l]][[k]]<-allTs[[va]] }    
}

pdf(paste("results/map_",labl,".pdf",sep=""),8,20)
for (l in (1:pm))  { 
  print(l)
  print(tmap_arrange(allT3[[l]],ncol=1,nrow=10)) }
dev.off()

save(allT3,allTs,file=paste("results/maptab_",labl,sep=""))

print(Sys.time())
print("successfully done")
q()



  

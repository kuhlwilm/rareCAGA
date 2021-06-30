#!/usr/bin/r

#module load xz TABIX/0.2.6 gcc/6.3.0 R/3.5.0 BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops; R --vanilla
library("GenomicRanges")
options("scipen"=100)
library(scales)
library(ggplot2)
require(gridExtra) 
library(rtracklayer)
library("gridGraphics")
'%ni%' <- Negate('%in%')
bps<-c("A","C","G","T")

## target space
trget<-read.table("/scratch/devel/cfontser/PanAf/CHR21/Final_Bam_August/Agilent_chr21_NR_chimptoHg19_sorted.bed",sep="\t",header=F)
trgets<-IRanges(start=trget[,2],end=trget[,3])

# sample names and metadata
load(file="/scratch/devel/mkuhlwilm/paf/snpad/colnams6.Robject")

# whole genome data was missing
snam4[((nrow(snam4)-58):nrow(snam4)),12]<-60
snam4[((nrow(snam4)-58):nrow(snam4)),21]<-"Chimp"
snam4[((nrow(snam4)-58):nrow(snam4)),13]<-"Keep"
snam4[((nrow(snam4)-58):nrow(snam4)),2]<-c(rep("EasternW",3),rep("Nigeria-CameroonW",10),rep("EasternW",15),rep("CentralW",19),rep("WesternW",12))
snam4[((nrow(snam4)-58):nrow(snam4)),4]<-"wGenom"
# select only relevant individuals: >1x coverage, <0.5% contamination, no PCA outlier, no excess heterozygosity, no uncertain subspecies origin
bind<-which(as.numeric(as.character(snam4[,12]))>1 & as.numeric(as.character(snam4[,17]))<.5 & snam4[,13]=="Keep"& !is.na(snam4[,1]) & snam4[,1]%ni%c("Kor1-35","CMNP1-24","Uga2-81","Gas1-10"))

# names of locations and samples
coln<-levels(as.factor(snam4[,2]))
sase<-snam4[bind,]
sasen<-snam4[bind,2]
coln<-coln[which(coln%in%sase[,2])]

save(coln,snam4,wind,file="/scratch/devel/mkuhlwilm/paf/snpad/coln")
wind=bind
print(length(bind))

yvec<-c(15:49);

# function to calculate frequency in each line
ffun<-function(typ,iput) {
  iput=iput[,which(sasen==typ),drop=F]
  apply(iput,1,frqfun)
}
frqfun<-function(input) {    if(length(which(is.na(input))==T)==length(input)) { return(NA) } else { return( length(which(input>0))/(length(!is.na(input))*2) ) }  }

# object to collect private sites
apri<-list();asites<-list()
#save(apri,asites,file="/scratch/devel/mkuhlwilm/paf/snpad/aprivat/aprivatT9")
#load(file="/scratch/devel/mkuhlwilm/paf/snpad/aprivat/aprivatT9")

for (yy in yvec) {
  print(yy)
  # load previously filtered medium-quality data
  tr<-try(load(file=paste("/scratch/devel/mkuhlwilm/paf/gt3_alldata_",yy,".Robject",sep="")))
  
  # get only target space SNPs
  nfo3<-    as.data.frame(subsetByOverlaps(IRanges(start=as.numeric(nfo2[,1]),end=as.numeric(nfo2[,1])),trgets))[,1]
  selec<-which(nfo2[,1]%in%nfo3)  
  nfo2<-nfo2[selec,]
  gtAQ<-gtAQ[selec,]
  gtAQ<-gtAQ[which(nfo2[,3]%in%bps),]
  nfo2<-nfo2[which(nfo2[,3]%in%bps),]
  
  ## only segregating sites within the reference panel
  rosu<-rowSums(matrix(as.numeric(gtAQ[,wind]),nrow=nrow(gtAQ[,wind])),na.rm=T)
  selsi<-which(rosu>0 )
  if(length(selsi)==0) { next }
  
  ipu<-gtAQ[selsi,wind]
  
  # get allele frequencies
  freqs<-lapply(coln,ffun,iput=ipu)
  freqs<-do.call(cbind,freqs)
  
  ## get almost private sites
  # somewhat permissive parameters: observed in one, cumulative frequency <1 in other sampling locations (take into account 2 chromosomes per individual)
  prifun<-function(typ,input) {   ifelse(input[,typ]>0 & rowSums(input[,-typ],na.rm=T)<0.5,1,ifelse(is.na(input[,typ]==T),NA,0)) }
  privat<-lapply(1:ncol(freqs),prifun,input=freqs)
  privat<-do.call(cbind,privat)
  apri[[yy]]<-privat
  asites[[yy]]<-nfo2[selsi,] ## include allele information for merging with other datasets
  # save the private alleles
  save(apri,asites,file="/scratch/devel/mkuhlwilm/paf/snpad/aprivat/aprivatT9")
}

print("done")
print(Sys.time())
q()


# write into bed file for easier access later on
load(file="/scratch/devel/mkuhlwilm/paf/snpad/aprivat/aprivatT9")
asitt<-do.call(rbind,asites)[,1]
for (yy in c(15:49)) { 
  obed<-cbind("chr21",as.numeric(asites[[yy]][,1])-1,as.numeric(asites[[yy]][,1]))
  write.table(obed,file=paste("/scratch/devel/mkuhlwilm/paf/snpad/aprivat/derisites/derived_9sites",yy,".bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
}



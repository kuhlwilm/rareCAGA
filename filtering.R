#!/usr/bin/r
yy=as.character(unlist((commandArgs(TRUE))))
print(yy)
print(Sys.time())

fifu<-function(input) { if(input!=".") { return(min(as.numeric(unlist(strsplit(input,split=","))))) } else { NA } }
bps<-c("A","C","G","T",".")
options("scipen"=100)

## this represents batches of 1 Mbp on chr21 (14M to 49M), to avoid memory failure
yy=as.numeric(yy)

yle=seq(0,49000000,1000000)
print (yy)
intvl<-c(yle[yy],yle[yy+1])

# use bcftools to get relevant information
agt<-system(paste("bcftools view -Ou -M2 -r chr21:",intvl[1],"-",intvl[2]," paf_chr21AF.vcf.gz | bcftools query -u -f '%POS %REF %ALT [%GT,] [%DP,] [%GQ,] [%AD-]\n' -",sep=""),intern=T)
if (length(agt)>0) {  
    agt<-do.call(rbind,strsplit(agt,split=" "))
    nfo<-agt[,1:3,drop=F]
    gt<-do.call(rbind,strsplit(agt[,4],split=","))
    dp<-matrix(as.numeric(do.call(rbind,strsplit(agt[,5],split=","))),nrow=nrow(nfo))
    gq<-matrix(as.numeric(do.call(rbind,strsplit(agt[,6],split=","))),nrow=nrow(nfo))
    ad<-matrix(do.call(rbind,strsplit(agt[,7],split="-")),nrow=nrow(nfo))
    ad<-matrix(as.numeric(unlist(apply(ad,c(1:2),fifu))),nrow=nrow(ad))
    
    # filtering for coverage 2-99 and genotype quality >=20, and, in case, allele balance
    GQ<-gt
    GQ[which(dp<2)]<-"./."
    GQ[which(dp>99)]<-"./."
    GQ[which(gq<20)]<-"./."
    GQ[which(gt%in%c("0/1","1/2","2/1","0/2","0/3","1/3","3/1","3/2","2/3","0/4","1/4","2/4","3/4","4/1","4/2","4/3") & (((ad/dp)<0.1)|ad==1))]<-"./."
    gtQ<-GQ    
    save(nfo,gtQ,file=paste("gt3_data_",yy,".Robject",sep=""))

    ## create genotype matrix
    # only use SNPs
    selc<-which(nfo[,2]%in%bps & nfo[,3]%in%bps)
    nfo2<-nfo[selc,]
    gtAQ<-gtQ[selc,]
    ## turn matrix into 0 1 2
    gtAQ[which(gtAQ=="1/1")]<-2;gtAQ[which(gtAQ=="0/1")]<-1;gtAQ[which(gtAQ=="0/0")]<-0;gtAQ[which(gtAQ=="./.")]<-NA
    save(gtAQ,nfo2,file=paste("gt3_alldata_",yy,".Robject",sep=""))
    }

q()



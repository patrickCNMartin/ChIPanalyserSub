###### RNA scaling ######

## Data set up starting with the usual

# Loading libraris and sourcing code
direc <- getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)


## sourcing scripts for analysis
## Loading source code files 
## Equivalent to:
library(ChIPanalyser)
setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)


## data
rna<-read.csv("../fpkm_dm3_cells.csv",stringsAsFactors=F) 

ctcf<-rna[grep("CTCF",rna[,6]),]
beaf<-rna[grep("BEAF",rna[,6]),]
suhw<-rna[grep("Hw",rna[,6]),]

ctcf<-ctcf[,which(grepl("BG3",colnames(ctcf))|grepl("Kc167",colnames(ctcf)) |grepl("S2",colnames(ctcf)))]

beaf<-beaf[,which(grepl("BG3",colnames(beaf))|grepl("Kc167",colnames(beaf)) |grepl("S2",colnames(beaf)))]

suhw<-suhw[,which(grepl("BG3",colnames(suhw))|grepl("Kc167",colnames(suhw)) |grepl("S2",colnames(suhw)))]



ctcf<-c(ctcf[2],ctcf[1])
beaf<-c(beaf[1],beaf[3])
suhw<-c(suhw[3],suhw[2])

scale<-list(ctcf,beaf,suhw)

### RNA re scale

 
 
 
### 

#  Using BEAF-32 modencode 922 S2 cell as template

# BEAF-32 BG3 modencode 921 as test 

## see above for scaling but ratio was obtain by tho folowing method 

# bm/S2mRNA * BG3mRNA

### RNA Scaling plots

#Data Loading 
AccessKc<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_Kc167.Rda"))
AccessBG3<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_BG3.Rda"))
AccessS2<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_S2.Rda"))

Access<-list(AccessBG3,AccessS2,AccessKc)
AccessOriginal<-list(AccessKc,AccessBG3,AccessS2)


DNASequenceSet<-get(load("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DNASequenceSet.Rda"))
input<-read.table("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputTRLPeaks.txt",sep=' ', comment.char='@',stringsAsFactors=F)
 
pfms<-list("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/CTCF.pfm","/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/BEAF-32.pfm","/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/su(Hw).pfm")

## rna shit
rna<-read.csv("../fpkm_dm3_cells.csv",stringsAsFactors=F) 

ctcf<-rna[grep("CTCF",rna[,6]),]
beaf<-rna[grep("BEAF",rna[,6]),]
suhw<-rna[grep("Hw",rna[,6]),]

ctcf<-ctcf[,which(grepl("BG3",colnames(ctcf))|grepl("Kc167",colnames(ctcf)) |grepl("S2",colnames(ctcf)))]

beaf<-beaf[,which(grepl("BG3",colnames(beaf))|grepl("Kc167",colnames(beaf)) |grepl("S2",colnames(beaf)))]

suhw<-suhw[,which(grepl("BG3",colnames(suhw))|grepl("Kc167",colnames(suhw)) |grepl("S2",colnames(suhw)))]



ctcf<-c(ctcf[2],ctcf[1])
beaf<-c(beaf[1],beaf[3])
suhw<-c(suhw[3],suhw[2])

scale<-list(ctcf,beaf,suhw)



method<-"MSEMean"
cores<-1
### Loading predicted vals 
## predicted in Kc cells and rescaled in BG3
ctcfKC<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_modEncode_908_CTCF_reduce50sigmoid_AUCMeanoptimalOutput.Rda"))
ctcfKCChIP<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_modEncode_908_CTCF_reduce50sigmoid_AUCMeanChIP.Rda"))

## predcited in BG3 and rescaled either in S2
beafBG3<-get(load("../../performAnalysis/peakWindowreduce20/BG3_modEncode_921_BEAF-32_reduce50sigmoid_AUCMeanoptimalOutput.Rda"))
beafBG3ChIP<-get(load("../../performAnalysis/peakWindowreduce20/BG3_modEncode_921_BEAF-32_reduce50sigmoid_AUCMeanChIP.Rda"))


## predicted in S2 and rescaled in Kc
suhwS2<-get(load("../../performAnalysis/peakWindowreduce20/S2_modEncode_331_Su(Hw)_reduce50sigmoid_AUCMeanoptimalOutput.Rda"))
suhwS2ChIP<-get(load("../../performAnalysis/peakWindowreduce20/S2_modEncode_331_Su(Hw)_reduce50sigmoid_AUCMeanChIP.Rda"))
 
 
predictions<-list(ctcfKC,beafBG3,suhwS2)
ChIPProfiles<-list(ctcfKCChIP,beafBG3ChIP,suhwS2ChIP)
 
 ### The one you will be working with or predicing in to show the resclaing 
 
ctcfBG3<-chipLoading("/home/pm16057/ChIP/modEncode_BG3/modEncode/modEncode_3674/signal_data_files/CTCF:Cell-Line=ML-DmBG3-c2#Developmental-Stage=Larvae-3rd-instar#RNAi-reagent=CG8573-RNAi#Tissue=CNS-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3674:repset.12058903.smoothedM.bed")

beafS2<-chipLoading("/home/pm16057/ChIP/modEncode_s2/modEncode/modEncode_922/signal_data_files/BEAF-32:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_922:repset.4620888.smoothedM.bed")

suhwKC<-chipLoading("/home/pm16057/ChIP/modEncode_s2/modEncode/modEncode_3719/signal_data_files/Su(Hw):Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#RNAi-reagent=Fly-LacZ-RNAi#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3719:repset.11436076.smoothedM.bed")
 
rescaledChIP<-list(ctcfBG3,beafS2,suhwKC) 

## subset of loci 

ctcfsub<-c("chr3L:6740001..6760000","chr3R:16760001..16780000")

beafsub<-c("chr3L:3800001..3820000","chrX:2160001..2180000")

suhwsub<-c("chr4:1120001..1140000","chr2R:19500001..19520000")

TFsub<-list(ctcfsub,beafsub,suhwsub)
 

## extract top parameters 

extract<-vector("list",3)
sets<-vector("list",3)

rescaledExtract<-vector("list",3)

topparam<-vector("list",3)

for(i in seq_along(predictions)){
    extract[[i]]<-vector("list", 4)
    rescaledExtract[[i]]<-vector("list",5)
    paraSet<-as.numeric(predictions[[i]][[1]][[1]][[method]])
    topparam[[i]]<-paraSet
    print(paraSet)
    extract[[i]][[1]]<-searchSites(predictions[[i]][[2]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
    extract[[i]][[2]]<-searchSites(predictions[[i]][[3]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
    extract[[i]][[3]]<-searchSites(predictions[[i]][[4]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
    extract[[i]][[4]]<-searchSites(predictions[[i]][[4]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
    
    ## sequence sets
    name<-names(ChIPProfiles[[i]])
 
    chr<-sapply(strsplit(name,":"),"[[",1)
    start<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",1)
    end<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",2)
    setSequence<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))
    names(setSequence)<-names(ChIPProfiles[[i]])
    sets[[i]]<-setSequence
    
   
    ChIPProfiles[[i]]<-ChIPProfiles[[i]][TFsub[[i]]]
    sets[[i]]<-setSequence[TFsub[[i]]]
    buff<-TFsub[[i]]
    #names(subSet)<-names(setSequence)[sub]
    extract[[i]][[1]]<-extract[[i]][[1]][[1]][buff]
    extract[[i]][[2]]<-extract[[i]][[2]][[1]][buff]
    extract[[i]][[3]]<-extract[[i]][[3]][[1]][buff]
  
    
 
    rescaledExtract[[i]][[1]]<-processingChIPseq(rescaledChIP[[i]],setSequence)
    OPP<-rescaledExtract[[i]][[1]][[2]]
    ChIPPrediction<-rescaledExtract[[i]][[1]][[1]][TFsub[[i]]]
    chips<-rescaledExtract[[i]][[1]][[1]]
    
                                                                        
    GPP<-genomicProfileParameters(PFM=pfms[[i]], PFMFormat="JASPAR", BPFrequency=DNASequenceSet,ScalingFactorPWM=paraSet[1])
    gw<-computeGenomeWidePWMScore(DNASequenceSet,GPP,Access[[i]])
    ## for extrcation sake 
    pwm<-computePWMScore(DNASequenceSet,gw,sets[[i]],Access[[i]])
    ###
    fullpwm<-computePWMScore(DNASequenceSet,gw,setSequence,Access[[i]])
    bm<-as.numeric(round((paraSet[2]/scale[[i]][[1]])*scale[[i]][[2]]))
    boundMolecules(OPP)<-bm
    rescaledExtract[[i]][[2]]<-computeOccupancy(pwm,Access[[i]],OPP)
    fulloccup<-computeOccupancy(fullpwm,Access[[i]],OPP)
    rescaledExtract[[i]][[3]]<-computeChipProfile(setSequence=sets[[i]],occupancyProfileParameters=OPP,occupancy=rescaledExtract[[i]][[2]])
    fullchip<-computeChipProfile(setSequence=setSequence,occupancyProfileParameters=OPP,occupancy=fulloccup)
    
    rescaledExtract[[i]][[4]]<-profileAccuracyEstimate(LocusProfile = ChIPPrediction, predictedProfile = rescaledExtract[[i]][[3]],occupancyProfileParameters = OPP,method="all")
    chips<-chips[!is.na(match(names(chips),names(fullchip[[1]])))]
    fullchip[[1]]<-fullchip[[1]][!is.na(match(names(fullchip[[1]]),names(chips)))]
    rescaledExtract[[i]][[5]]<-profileAccuracyEstimate(LocusProfile = chips, predictedProfile = fullchip,occupancyProfileParameters = OPP,method="all")
    
    


}



noaccessS2<-.AccessExtract(subSet,AccessS2)
noaccessBG3<-.AccessExtract(subSet,AccessBG3)
noaccessKc167<-.AccessExtract(subSet,AccessKc)
 
tf<-c("CTCF","BEAF-32","Su(Hw)")
cell<-c("Kc167","BG3","S2","BG3","S2","Kc167")
labs<-LETTERS[1:6]

pdf(paste0("RNA_rescale_ChIPanalyser_","MSE",".pdf"), width=20,height=15)
#par(oma=c(0,0,9,0))
layout(matrix(cbind(c(1,2,5,6,9,10),c(3,4,7,8,11,12),c(13,13,14,14,15,15)),ncol=3), width=c(7,7,3),height=c(1,1,1))
par(family="mono")
par(xpd=T)

cols<-c("#4f9da6","#233142","#ff5959","#facf5a")
count<-0
for(i in seq_along(rescaledExtract)){
    ## orginal data 
    count<-count+1
    for(k in seq_along(sets[[i]])){
       x<-seq(start(sets[[i]])[k],end(sets[[i]])[k],by=10)
       x<-c(x[1]-1,x,x[length(x)]+1)
       if(k==1){par(mar=c(4,2,4.5,2))}else{par(mar=c(4,2,3.8,2))}
    
       plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(sets[[i]])[k],end(sets[[i]])[k]),ylim=c(0,1))
       
       title(xlab=paste0("Genomic Position on ",as.character(seqnames(sets[[i]]))[k]),cex.lab=1.5)
       
       if(k==1){
       title(main=paste0(tf[i]," in ",cell[i]," - lambda = ",topparam[[i]][1]," & Bound Molecules = ",topparam[[i]][2]),cex.main=1.8,line=0.5)
       #text(x=(x-1000), y=1.2, labels=labs[count])
       }
       
       axis(1,at=round(seq(start(sets[[i]])[k],end(sets[[i]])[k],length.out=10)),labels=round(seq(start(sets[[i]])[k],end(sets[[i]])[k],length.out=10)),cex.axis=1.5)
       
       noaccess<-.AccessExtract(sets[[i]],AccessOriginal[[i]])
       
       for(j in seq_len(nrow(noaccess[[k]]))){
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=45,lwd=1,border=NA)
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
       }
       
       
       chipInd<-c(0,ChIPProfiles[[i]][[k]][seq(0,length(ChIPProfiles[[i]][[k]]),by=10)],0)
       predInd<-c(0,extract[[i]][[2]][[k]]$ChIP,0)

       polygon(x,chipInd,density=NA,col="#233142",lwd=2)
       lines(x,predInd,col="#ff5959",lwd=2.5)
       occupancy<-extract[[i]][[1]][[k]]

      OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), round(0.9*length(occupancy$Occupancy)))]

      ReScale<-((OccupScaling$Occupancy/max(OccupScaling$Occupancy)))
  
      lines(x=start(OccupScaling),y=ReScale*0.8,type="h",col="#4f9da6",lwd=1.8)

    
     }
     
     # rescaled data 
     
     for(k in seq_along(sets[[i]])){
       x<-seq(start(sets[[i]])[k],end(sets[[i]])[k],by=10)
       x<-c(x[1]-1,x,x[length(x)]+1)
        if(k==1){par(mar=c(4,2,4.5,2))}else{par(mar=c(4,2,3.8,2))}
    
       plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(sets[[i]])[k],end(sets[[i]])[k]),ylim=c(0,1))
       title(xlab=paste0("Genomic Position on ",as.character(seqnames(sets[[i]]))[k]),cex.lab=1.5)
       
       axis(1,at=round(seq(start(sets[[i]])[k],end(sets[[i]])[k],length.out=10)),labels=round(seq(start(sets[[i]])[k],end(sets[[i]])[k],length.out=10)),cex.axis=1.5)
       if(k==1){
       title(main=paste0(tf[i]," in ", cell[i+3]," - ",names(rescaledExtract[[i]][[3]])),cex.main=1.8,line=0.5)
       #text(x=(x-1000), y=1.2, labels=labs[count])
       }
      
       noaccess<-.AccessExtract(sets[[i]],Access[[i]])
       
       for(j in seq_len(nrow(noaccess[[k]]))){
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=45,lwd=1,border=NA)
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
           rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
       }
       
       bufferChIP<-rescaledExtract[[i]][[1]][[1]][TFsub[[i]]]
       chipInd<-c(0,bufferChIP[[k]][seq(0,length(bufferChIP[[k]]),by=10)],0)
       
       predBuffer<-rescaledExtract[[i]][[3]][[1]][[k]]
       predInd<-c(0,predBuffer$ChIP,0)

       polygon(x,chipInd,density=NA,col="#233142",lwd=2)
       lines(x,predInd,col="#ff5959",lwd=2.5)
       occupancy<- AllSitesAboveThreshold(rescaledExtract[[i]][[2]])[[1]][[k]]

      OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), round(0.9*length(occupancy$Occupancy)))]

      ReScale<-((OccupScaling$Occupancy/max(OccupScaling$Occupancy)))
  
      lines(x=start(OccupScaling),y=ReScale*0.8,type="h",col="#4f9da6",lwd=1.8)

    
     }
    
   

}






## initating boxplot for loop
labs<-LETTERS[7:9]
for(i in seq_along(rescaledExtract)){
    ext<-extract[[i]][[4]][[1]]
    meth<-gsub("Mean","",method)
    ext<-sapply(ext,function(x){x[[1]][meth]})
    rescale<-rescaledExtract[[i]][[5]][[1]]
    rescale<-sapply(rescale,function(x){x[[1]][meth]})
    
    dat<-list("Estimated"=ext,"Rescaled"=rescale)
    boxplot(dat,main="MSE Distribution",col=c("#4f9da6","#233142"),frame=F,cex.axis=1.5,cex.main=1.8)
    text(x=(-0.1),y=0.055, labels=labs[i]) 
}


dev.off()  
    
  

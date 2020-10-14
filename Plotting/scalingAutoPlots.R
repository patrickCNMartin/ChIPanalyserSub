
direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm3)
library(BSgenome)
library(RcppRoll)
library(GenomicRanges)
library(ROCR)
library(zoo)


setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyser_1.1")
files <- dir()
for (i in files) source(i)
setwd(direc)








 
 




name<-c("recallMean","FscoreMean","MCCMean","AUCMean")

for( method in name){

ChIPEstimation<-get(load("Kc167_modEncode_908_CTCF_reduce20sigmoidChIP.Rda"))

load("Kc167_modEncode_908_CTCF_reduce20sigmoidoptimalOutput.Rda")

AccessKc<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_Kc167.Rda"))
 AccessBG3<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_BG3.Rda"))
  AccessS2<-get(load("~/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_S2.Rda"))
  DNASequenceSet<-get(load("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DNASequenceSet.Rda"))
  input<-read.table("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputTRLPeaks.txt",sep=' ', comment.char='@',stringsAsFactors=F)
 
 ### The one you will be working with or predicing in 
 ChIPPrediction<-chipLoading("/home/pm16057/ChIP/modEncode_BG3/modEncode/modEncode_3674/signal_data_files/CTCF:Cell-Line=ML-DmBG3-c2#Developmental-Stage=Larvae-3rd-instar#RNAi-reagent=CG8573-RNAi#Tissue=CNS-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3674:repset.12058903.smoothedM.bed")
 pfm<-"/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/CTCF.pfm"
optimalParama <- get(load("optimal.Rda"))
  
  
paraSet<- as.numeric(optimalParama[[1]][[method]])

     occup <- searchSites(optimal[[2]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
     chip<-searchSites(optimal[[3]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
     preds<-searchSites(optimal[[4]],ScalingFactor=paraSet[1],BoundMolecules=paraSet[2])
 occupEstimation<-occup
 chipEstimation<-chip
 predsEstimation<-preds
 name<-names(ChIPEstimation)
 
 chr<-sapply(strsplit(name,":"),"[[",1)
 start<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",1)
 end<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",2)
 setSequence<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))
 names(setSequence)<-names(ChIPEstimation)
 
 
 setSequence<-setSequence[names(setSequence) %in% names(chip[[1]])] 
  ChIPEstimation<-ChIPEstimation[names(ChIPEstimation) %in% names(chip[[1]])] 
 #### Find the ones you want to keep
 sub<-c(4,8,14,16)
 
 
 ChIPEstimation<-ChIPEstimation[sub]
 subSet<-setSequence[sub]
#names(subSet)<-names(setSequence)[sub]
occupEstimation<-occupEstimation[[1]][sub]
chipEstimation<-chipEstimation[[1]][sub]
 predsEstimation<-predsEstimation[[1]][sub]

 
 ChIPPrediction<-processingChIPseq(ChIPPrediction,subSet)
OPP<-ChIPPrediction[[2]]
ChIPPrediction<-ChIPPrediction[[1]]

 


 


GPP<-genomicProfileParameters(PFM=pfm, PFMFormat="JASPAR", BPFrequency=DNASequenceSet,ScalingFactorPWM=paraSet[1])
gw<-computeGenomeWidePWMScore(DNASequenceSet,GPP,AccessBG3)
pwm<-computePWMScore(DNASequenceSet,gw,subSet,AccessBG3)

##rescaling 
bm<-as.numeric(round((paraSet[2]/beaf[3])*beaf[1]))
boundMolecules(OPP)<-bm
occupPrediction<-computeOccupancy(pwm,AccessBG3,OPP)
names(subSet)<-names(AllSitesAboveThreshold(occupPrediction)[[1]])
 chipPrediction<-computeChipProfile(setSequence=subSet,occupancyProfileParameters=OPP,occupancy=occupPrediction)[[1]]
occupPrediction<-AllSitesAboveThreshold(occupPrediction)
occupPrediction<-occupPrediction[[1]]


noaccessS2<-.AccessExtract(subSet,AccessS2)
 noaccessBG3<-.AccessExtract(subSet,AccessBG3)
 noaccessKc167<-.AccessExtract(subSet,AccessKc)

 

## you need to set up your data, which is not here because you are an idiot :D
pdf(paste0("reScale_CTCF_Kc167_to_BG3",method,".pdf"), width=20,height=10)
par(oma=c(0,0,9,0))
par(mfrow=c(4,2))
par(family="mono")
par(xpd=T)
#cols<-c("#6b7d81","#475263","#970e47","#051727")
#cols<-c("#6b7d81","#475263","#a4ced4","#051727")
cols<-c("#F0E442","#999999","#D55E00","#56B4E9")

for(i in seq_along(ChIPPrediction)){
    
    x<-seq(start(subSet)[i],end(subSet)[i],by=10)
    x<-c(x[1]-1,x,x[length(x)]+1)
    
    ### estimated 
    par(mar=c(4,2,2,2))
    
    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(subSet)[i],end(subSet)[i]),ylim=c(0,1))
    title(xlab=paste0("Genomic Position on ",as.character(seqnames(subSet))[i]),cex.lab=1.3)
   
    if(i ==1){
    title(main=paste0("Estimated Parameters in Kc167 Cells: lambda = ",paraSet[1]," & N = ",paraSet[2]),cex.main=1.6) 
    text(start(subSet)[i],1.1,labels="A", font=2,cex=3)
    }
    
    axis(1,at=round(seq(start(subSet)[i],end(subSet)[i],length.out=10)),labels=round(seq(start(subSet)[i],end(subSet)[i],length.out=10)),cex.axis=1.2)
    
   
    
     for(j in seq_len(nrow(noaccessKc167[[i]]))){
    rect(noaccessKc167[[i]][j,"start"],0,noaccessKc167[[i]][j,"end"],0.9,col=cols[1],density=10,angle=45,lwd=0.8,border=NA)
    rect(noaccessKc167[[i]][j,"start"],0,noaccessKc167[[i]][j,"end"],0.9,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
    rect(noaccessKc167[[i]][j,"start"],0,noaccessKc167[[i]][j,"end"],0.9,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
    }


    chipInd<-c(0,ChIPEstimation[[i]][seq(0,length(ChIPEstimation[[i]]),by=10)],0)
    predInd<-c(0,chipEstimation[[i]]$ChIP,0)

    polygon(x,chipInd,density=NA,col=cols[2],lwd=2)
    lines(x,predInd,col=cols[3],lwd=2)
    occupancy<-occupEstimation[[i]]

    OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), round(0.9*length(occupancy$Occupancy)))]

    ReScale<-((OccupScaling$Occupancy/max(OccupScaling$Occupancy)))
    #sts<-start(OccupScaling)
    #for( k in seq_along(sts)){
        #lines(x=c(sts[k],sts[k]),y=c(0,ReScale[k]),lwd=2, col=cols[4])
   # }
    lines(x=start(OccupScaling),y=ReScale*0.8,type="h",col=cols[4],lwd=1.5)
    
    
    ### Rescale
    
    par(mar=c(4,2,2,4))
    
    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(subSet)[i],end(subSet)[i]),ylim=c(0,1))
    title(xlab=paste0("Genomic Position on ",as.character(seqnames(subSet))[i]),cex.lab=1.3)
   
    if(i ==1){title(main=paste0("Fitted Parameters in BG3 Cells: lambda = ",paraSet[1]," & N = ",bm),cex.main=1.6) 
       text(start(subSet)[i],1.1,labels="B", font=2,cex=3)
    }
    
    axis(1,at=round(seq(start(subSet)[i],end(subSet)[i],length.out=10)),labels=round(seq(start(subSet)[i],end(subSet)[i],length.out=10)),cex.axis=1.2)
    
    
     for(j in seq_len(nrow(noaccessBG3[[i]]))){
    rect(noaccessBG3[[i]][j,"start"],0,noaccessBG3[[i]][j,"end"],0.9,col=cols[1],density=10,angle=45,lwd=0.8,border=NA)
    rect(noaccessBG3[[i]][j,"start"],0,noaccessBG3[[i]][j,"end"],0.9,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
    rect(noaccessBG3[[i]][j,"start"],0,noaccessBG3[[i]][j,"end"],0.9,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
    }

    chipInd<-c(0,ChIPPrediction[[i]][seq(0,length(ChIPPrediction[[i]]),by=10)],0)
    predInd<-c(0,chipPrediction[[i]]$ChIP,0)

    polygon(x,chipInd,density=NA,col=cols[2],lwd=2)
    lines(x,predInd,col=cols[3],lwd=2)
    occupancy<-occupPrediction[[i]]

    OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), round(0.9*length(occupancy$Occupancy)))]

    ReScale<-((OccupScaling$Occupancy/max(OccupScaling$Occupancy)))
     lines(x=start(OccupScaling),y=ReScale,type="h",col=cols[4],lwd=1.5)
    
    

}
mtext("Cell Line Prediction for CTCF",outer=T,cex=2.2,line=4)

dev.off()


}




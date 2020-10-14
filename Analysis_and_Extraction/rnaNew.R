###### RNA scaling ######


### RNA re scale

direc <- "/home/pm16057/ChIPanalyser/ChIPanalyserFinal/RNA/refine"

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)


## sourcing scripts for analysis
setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)



###

#  Using BEAF-32 modencode 922 S2 cell as template
# optimal Paramters were lambda =5 and BM= 50000 / DHS access 500bp /GeometricMean
# BEAF-32 BG3 modencode 921 as test

## see above for scaling but ratio was obtain by tho folowing method

# bm/S2mRNA * BG3mRNA

### RNA Scaling plots

#Data Loading
AccessKc<-get(load("/home/pm16057/DNAaccess/cellAccess/Kc_DHS_005.Rda"))
AccessBG3<-get(load("/home/pm16057/DNAaccess/cellAccess/BG3_DHS_005.Rda"))
AccessS2<-get(load("/home/pm16057/DNAaccess/cellAccess/S2_DHS_005.Rda"))

Access<-list("BG3"=AccessBG3,"S2"=AccessS2,"S2"=AccessS2)
AccessOriginal<-list(AccessKc,AccessBG3,AccessKc)


DNASequenceSet<-get(load("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DNASequenceSet.Rda"))
#input<-read.table("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputTRL.txt",sep=' ', comment.char='@',stringsAsFactors=F)

pfms<-list("CTCF"="/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/CTCF.pfm",
           "beaf"="/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/BEAF-32.pfm",
           "hw"="/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/su(Hw).pfm")

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
suhw<-c(suhw[2],suhw[3])

scale<-list("CTCF"=ctcf,"beaf"=beaf,"hw"=suhw)



method<-"MSE"
cores<-1
### Loading predicted vals
## predicted in Kc cells and rescaled in BG3
ctcfKC<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_GSM762842_CTCF_reduce_top10opti3293_Kc_DHS_005sigmoid_OptimalOutputTraining.Rda"))
ctcfKCChIP<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_GSM762842_CTCF_reduce_top10opti3293_Kc_DHS_005sigmoid_ChIPTraining.Rda"))

## predcited in BG3 and rescaled either in S2
beafBG3<-get(load("../../performAnalysis/peakWindowreduce20/BG3_modEncode_921_BEAF-32_reduce_top10opti3293_BG3_DHS_005sigmoid_OptimalOutputTraining.Rda"))
beafBG3ChIP<-get(load("../../performAnalysis/peakWindowreduce20/BG3_modEncode_921_BEAF-32_reduce_top10opti3293_BG3_DHS_005sigmoid_ChIPTraining.Rda"))


## predicted in S2 and rescaled in Kc
suhwKc<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_Su(Hw)_reduce_top10opti3293_Kc_DHS_005sigmoid_OptimalOutputTraining.Rda"))
suhwKcChIP<-get(load("../../performAnalysis/peakWindowreduce20/Kc167_Su(Hw)_reduce_top10opti3293_Kc_DHS_005sigmoid_ChIPTraining.Rda"))


predictions<-list("ctcf_Kc167"=ctcfKC,"beaf_BG3"=beafBG3,"hw_Kc167"=suhwKc)
ChIPProfiles<-list("ctcf_Kc167"=ctcfKCChIP,"beaf_BG3"=beafBG3ChIP,"hw_Kc167"=suhwKcChIP)

 ### The one you will be working with or predicing in to show the resclaing

ctcfBG3<-get(load("/home/pm16057/ChIP/modEncode_BG3/modEncode/modEncode_3674/signal_data_files/CTCF:Cell-Line=ML-DmBG3-c2#Developmental-Stage=Larvae-3rd-instar#RNAi-reagent=CG8573-RNAi#Tissue=CNS-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3674:repset.12058903.smoothedM.bed.Rda"))

beafS2<-get(load("/home/pm16057/ChIP/modEncode_s2/modEncode/modEncode_922/signal_data_files/BEAF-32:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_922:repset.4620888.smoothedM.bed.Rda"))

suhwS2<-get(load("/home/pm16057/ChIP/modEncode_s2/modEncode/modEncode_331/signal_data_files/Su(Hw):Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_331:repset.4621850.smoothedM.bed.Rda"))

rescaledChIP<-list("ctcf_BG3"=ctcfBG3,"beaf_S2"=beafS2,"hw_S2"=suhwS2)

## subset of loci

ctcfsub<-c(2,6)
beafsub<-c(1,6)

suhwsub<-c(6,8)

TFsub<-list("CTCF"=ctcfsub,"beaf"=beafsub,"hw"=suhwsub)

cores<-10
method<-"MSE"

rescaledExtract <- vector("list",3)
names(rescaledExtract)<-paste0(names(predictions),"_",names(rescaledChIP))

rescaledChIPsignal <- vector("list",3)
names(rescaledExtract)<-paste0(names(predictions),"_",names(rescaledChIP))

for(i in seq_along(predictions)){
      print(names(rescaledExtract)[i])
      ## first lets load data
      localOptimal<-predictions[[i]]
      localChIP <-ChIPProfiles[[i]]

      ## check point
      print(names(predictions)[[i]]==names(ChIPProfiles)[[i]])

      ## extract top paramters
      param <- as.numeric(localOptimal[[1]][[1]][[method]])
      print(param)

      lambda <- param[1]
      bound <- param[2]

      ## extracting and pre processing chip
      setSequence <- loci(localChIP)

      localExtractedChIP <- processingChIP(rescaledChIP[[i]],setSequence, cores=cores)
      rescaledChIPsignal[[i]]<-localExtractedChIP
      ## sclaing bound Molecules


      bmSet <- c(bound , as.numeric(round((bound/scale[[i]][[1]])*scale[[i]][[2]])),
                 bound/100)
      print(bmSet)
      rescaling<-list("Carry"=NULL,"Rescaled"=NULL,"100Fold"=NULL)
      for(j in seq_along(bmSet)){
          print(names(rescaling)[j])
          GPP <- genomicProfiles(PFM=pfms[[i]], PFMFormat="JASPAR", BPFrequency=DNASequenceSet,lambdaPWM=lambda,boundMolecules=bmSet[j],stepSize=100)
          gw <- computeGenomeWideScores(GPP,DNASequenceSet,Access[[i]],cores=cores)
          pwm <- computePWMScore(gw,DNASequenceSet,localExtractedChIP,Access[[i]],cores=cores)
          occup <- computeOccupancy(pwm)
          chip <- computeChIPProfile(occup,localExtractedChIP, cores=cores)
          GoF <- profileAccuracyEstimate(chip,localExtractedChIP,cores=cores)
          rescaling[[j]]<-list("Occupancy"=occup,"ChIPPrediction"=chip,"GoF"=GoF)

      }
      rescaledExtract[[i]]<-rescaling

}



## let's plots these mofos

## lets test dem profiles
#pdf("plotEverything.pdf", height=7, width =14)
#par(mar=c(2,2,2,2),mfrow=c(5,2))
#for(i in seq_along(rescaledExtract)){
#   for(j in seq_along(rescaledExtract[[i]])){
#       plotOccupancyProfile(rescaledExtract[[i]][[j]]$ChIPPrediction,ChIPProfiles[[i]],chromatinState=Access[[i]])
#   }
#}

#dev.off()

##### This works fine so that's neat

###

### Lets do the actual plotting
#tf<-c("CTCF","BEAF-32","su(Hw)")
#cell<-c("Kc167","BG3","Kc167","BG3","S2","S2")
#labs<-LETTERS[1:6]
## you need to set up your data, which is not here because you are an idiot :D
#pdf(paste0("RNA_rescale_ChIPanalyser_","MSE",".pdf"), width=27,height=15)
#par(oma=c(0,0,9,0))
#layout(matrix(cbind(c(1,2,5,6,9,10),c(3,4,7,8,11,12),c(13,13,14,14,15,15)),ncol=3), width=c(7,7,5.5),height=c(1,1,1))
#par(family="sans")
#par(xpd=NA)

#cols<-c("#ff5959","#233142","#facf5a","#facf5a")
#count<-0

#for(i in seq_along(rescaledExtract)){
#    setTrain <- loci(ChIPProfiles[[i]])
#    scoresTrain<- scores(ChIPProfiles[[i]])
#
#    param <- as.numeric(predictions[[i]][[1]][[1]][[method]])
#
#    lambda <- param[1]
#    bound <- param[2]
#    count<-count+1
#
#  for(k in TFsub[[i]]){
#    print(k)
#    predictionTrain <- searchSites(predictions[[i]]$ChIPProfiles, lambda, bound,names(scoresTrain)[k])
#    x<-seq(start(setTrain)[k],end(setTrain)[k],by=100)
#    x<-c(x[1]-1,x,x[length(x)]+1)
#    if(k==1){par(mar=c(4,2,4.5,2))}else{par(mar=c(4,2,3.8,2))}
#
#    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(setTrain)[k],end(setTrain)[k]),ylim=c(0,1))
#
#    title(xlab=paste0("Genomic Position on ",as.character(seqnames(setTrain))[k]),cex.lab=1.5)
#
#    if(k==1){
#    title(main=paste0(tf[i]," in ",cell[i]," - lambda = ",lambda," & Bound Molecules = ",bound),cex.main=1.8,line=0.5)
#    text(x=(x[1]-1000), y=1.2, labels=labs[count],cex=4)
#    }
#
#    axis(1,at=round(seq(start(setTrain)[k],end(setTrain)[k],length.out=10)),labels=round(seq(start(setTrain)[k],end(setTrain)[k],length.out=10)),cex.axis=1.5)
#
#    noaccess<-.AccessExtract(setTrain[k],AccessOriginal[[i]])[[1]]
#
#    for(j in seq_len(nrow(noaccess))){
#        rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
#        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
#        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
#    }
#
#     local<-scoresTrain[[k]]
#    chipInd<-c(0,local[seq(0,length(local),by=100)],0)
#    predInd<-c(0,predictionTrain[[1]][[1]]$ChIP,0)
#
#    polygon(x,chipInd,density=NA,col="#233142",lwd=2)
#    lines(x,predInd,col="#ff5959",lwd=2.5)
#  }
#  for(k in TFsub[[i]]){
#    print(k)
#    validationScore <- scores(rescaledChIPsignal[[i]])
#    validationSet <-loci(rescaledChIPsignal[[i]])
#
#    x<-seq(start(validationSet)[k],end(validationSet)[k],by=100)
#    x<-c(x[1]-1,x,x[length(x)]+1)
#    if(k==1){par(mar=c(4,2,4.5,2))}else{par(mar=c(4,2,3.8,2))}
#
#    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(validationSet)[k],end(validationSet)[k]),ylim=c(0,1))
#
#    title(xlab=paste0("Genomic Position on ",as.character(seqnames(validationSet))[k]),cex.lab=1.5)
#
#    if(k==1){
#    title(main=paste0(tf[i]," in ",cell[i]," - lambda = ",lambda," & Bound Molecules = ",bound),cex.main=1.8,line=0.5)
#    text(x=(x[1]-1000), y=1.2, labels=labs[count],cex=4)
#    }
#
#    axis(1,at=round(seq(start(validationSet)[k],end(validationSet)[k],length.out=10)),labels=round(seq(start(validationSet)[k],end(validationSet)[k],length.out=10)),cex.axis=1.5)
#
#    noaccess<-.AccessExtract(validationSet[k],AccessOriginal[[i]])[[1]]
#
#    for(j in seq_len(nrow(noaccess))){
#        rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
#        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
#        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
#    }
#
#     local<-validationScore[[k]]
#    chipInd<-c(0,local[seq(0,length(local),by=100)],0)
#    Carry<-profiles(rescaledExtract[[i]][[1]]$ChIPPrediction)
#    CarryInd<-c(0,Carry[[1]][[k]]$ChIP,0)
#
#    rescale <-profiles(rescaledExtract[[i]][[2]]$ChIPPrediction)
#    rescaleInd<-c(0,rescale[[1]][[k]]$ChIP,0)
#
#    fold <-profiles(rescaledExtract[[i]][[3]]$ChIPPrediction)
#    foldInd<-c(0,fold[[1]][[k]]$ChIP,0)
#
#
#    polygon(x,chipInd,density=NA,col="#233142",lwd=2)
#    lines(x,CarryInd,col="#56B4E9",lwd=2.5)
#    lines(x,rescaleInd,col="#ff5959",lwd=2.5,lty=2)
#    lines(x,foldInd,col="#009E73",lwd=2.5,lty=1)
#    legend(x=x[length(x)],y=0.75,col=c("#56B4E9","#ff5959","#009E73"),
#    lty=c(1,2,1),legend=c("Carry-Over","Rescaled","100 Fold"),bty="n",lwd=rep(2.5,3),cex=1.2)
#  }
#
#
#}#
#





## initating boxplot for loop
#labs<-LETTERS[7:9]
#par(xpd=NA)
#for(i in seq_along(rescaledExtract)){
#    ext<-extract[[i]][[4]][[1]]
#    meth<-gsub("Mean","",method)
#    ext<-sapply(ext,function(x){x[meth]})
#    rescale<-profiles(rescaledExtract[[i]][[5]])[[1]]
#    rescale<-sapply(rescale,function(x){x[meth]})
#    NoScale<-profiles(rescaledExtract[[i]][[8]])[[1]]
#    NoScale<-sapply(NoScale,function(x){x[meth]})
#    NoScale1000<-profiles(rescaledExtract[[i]][[11]])[[1]]
#    NoScale1000<-sapply(NoScale1000,function(x){x[meth]})
#    par(xpd=NA)
#    par(mar=c(4,7,4,2))
#    dat<-list("Estimated"=ext,"Rescaled"=rescale,"Carry-Over"=NoScale,"Carry-Over/1000"=NoScale1000)
#    boxplot(dat,main="MSE Distribution",col=c("#4f9da6","#ff5959","#facf5a","#009E73"),frame=F,cex.axis=1.5,cex.main=1.8,ylim=c(0,0.1))
#    text(x=(-0.25),y=0.11, labels=labs[i],cex=4)
#}


#dev.off()

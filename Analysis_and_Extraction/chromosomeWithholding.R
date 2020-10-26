###############################################################################
########################### Method comp #######################################
###############################################################################

direc <- getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
library(ChIPanalyser)

## sourcing scripts for analysis



## going to select regions from different chromosome
## we can base this on method data sets


#input <-read.table("~/ChIPanalyser/ChIPanalyserFinal/methods/DataInputMethod.txt", sep=' ', comment.char='@', stringsAsFactors=F)
input<-read.table("~/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputDHSgeometric.txt", sep=' ', comment.char='@', stringsAsFactors=F)


DNASequenceSet <- get(load(input[1,6]))
AccessBEAF<-get(load(input[7,7]))
AccessCTCF<-get(load(input[7,7]))
Accesssuhw<-get(load(input[7,7]))

## test data set clean shifted to beaf-32 instead just becase its cleaner

fileBEAF <- as.vector(as.matrix(input[7,5]))
fileCTCF <- as.vector(as.matrix(input[12,5]))
filesuhw <- as.vector(as.matrix(input[10,5]))
## loading loci

locus <- get(load(as.vector(as.matrix(input[7,8]))))

# laoding pfms

pfmBEAF<-as.vector(as.matrix(input[7,3]))
pfmCTCF<-as.vector(as.matrix(input[12,3]))
pfmsuhw<-as.vector(as.matrix(input[10,3]))

PFMFor <-"JASPAR"
## filter chromosomes

cores<-10

TrainingLoci <- locus[which(as.character(seqnames(locus)) %in% c("chr3R"))]
ValidationLoci <- locus[which(as.character(seqnames(locus)) %in% c("chr2R"))]
ValidationArtLoci <-locus[which(as.character(seqnames(locus)) %in% c("chr3R"))]

PO<-parameterOptions(noiseFilter="sigmoid",stepSize=100)




chipTrainBEAF <- processingChIP(fileBEAF,TrainingLoci, reduce=10,parameterOptions=PO,cores=cores)

chipValidationBEAF <- processingChIP(fileBEAF,ValidationLoci, reduce=20,
                                parameterOptions=PO, cores=cores)
chipValidationArtBEAF <- processingChIP(fileBEAF,ValidationArtLoci, reduce=20,
                                      parameterOptions=PO, cores=cores)

.scores(chipValidationArtBEAF)<-scores(chipValidationArtBEAF)[11:20]
.loci(chipValidationArtBEAF)<-loci(chipValidationArtBEAF)[11:20]




chipTrainCTCF <- processingChIP(fileCTCF,TrainingLoci, reduce=10,parameterOptions=PO,cores=cores)

chipValidationCTCF <- processingChIP(fileCTCF,ValidationLoci, reduce=20,
                                parameterOptions=PO, cores=cores)
chipTrainsuhw <- processingChIP(filesuhw,TrainingLoci, reduce=10,parameterOptions=PO,cores=cores)
chipValidationsuhw <- processingChIP(filesuhw,ValidationLoci, reduce=20,parameterOptions=PO, cores=cores)



GPBEAF <- genomicProfiles(PFM=pfmBEAF,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,stepSize=100)
GPCTCF <- genomicProfiles(PFM=pfmCTCF,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,stepSize=100)
GPsuhw <- genomicProfiles(PFM=pfmsuhw,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,stepSize=100)




## training

optimalBEAF <- computeOptimal(GPBEAF,DNASequenceSet,chipTrainBEAF,AccessBEAF,PO, cores=cores)
optimalCTCF <- computeOptimal(GPCTCF,DNASequenceSet,chipTrainCTCF,AccessCTCF,PO, cores=cores)
optimalsuhw <- computeOptimal(GPsuhw,DNASequenceSet,chipTrainsuhw,Accesssuhw,PO, cores=cores)

paramBEAF<-optimalBEAF[[1]][[1]][["MSE"]]
paramCTCF<-optimalCTCF[[1]][[1]][["MSE"]]
paramsuhw<-optimalsuhw[[1]][[1]][["MSE"]]

subOptiBEAF <- searchSites(optimalBEAF, paramBEAF[1],paramBEAF[2])
subOptiCTCF <- searchSites(optimalCTCF, paramCTCF[1],paramCTCF[2])
subOptisuhw <- searchSites(optimalsuhw, paramsuhw[1],paramsuhw[2])

pdf("BEAFtest.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(subOptiBEAF$ChIPProfiles, chipTrainBEAF, chromatinState=AccessBEAF)
dev.off()

pdf("CTCFtest.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(subOptiCTCF$ChIPProfiles, chipTrainCTCF, chromatinState=AccessCTCF)
dev.off()

pdf("suhwtest.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(subOptisuhw$ChIPProfiles, chipTrainsuhw, chromatinState=Accesssuhw)
dev.off()

## Validation

GPvalBEAF <-genomicProfiles(PFM=pfmBEAF,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,
                        stepSize=100,lambdaPWM=paramBEAF[1],boundMolecules=paramBEAF[2],stepSize=100)

GPvalCTCF <-genomicProfiles(PFM=pfmCTCF,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,
                        stepSize=100,lambdaPWM=paramCTCF[1],boundMolecules=paramCTCF[2])

GPvalsuhw <-genomicProfiles(PFM=pfmsuhw,PFMFormat=PFMFor, BPFrequency=DNASequenceSet,
                        stepSize=100,lambdaPWM=paramsuhw[1],boundMolecules=paramsuhw[2])

gwBEAF<-computeGenomeWideScores(GPvalBEAF,DNASequenceSet,AccessBEAF,cores=cores)
gwCTCF<-computeGenomeWideScores(GPvalCTCF,DNASequenceSet,AccessCTCF,cores=cores)
gwsuhe<-computeGenomeWideScores(GPvalsuhw,DNASequenceSet,Accesssuhw,cores=cores)

pwmBEAF<-computePWMScore(gwBEAF,DNASequenceSet,chipValidationBEAF,AccessBEAF,cores=cores)
pwmBEAFArt<-computePWMScore(gwBEAF,DNASequenceSet,chipValidationArtBEAF,AccessBEAF,cores=cores)


pwmCTCF<-computePWMScore(gwCTCF,DNASequenceSet,chipValidationCTCF,Access,cores=cores)
pwmsuhw<-computePWMScore(gwsuhe,DNASequenceSet,chipValidationsuhw,Access,cores=cores)

occupBEAF<-computeOccupancy(pwmBEAF)
occupBEAFArt<-computeOccupancy(pwmBEAFArt)


occupCTCF<-computeOccupancy(pwmCTCF)
occupsuhw<-computeOccupancy(pwmsuhw)

chipBEAF<-computeChIPProfile(occupBEAF,chipValidationBEAF,cores=cores)
chipBEAFArt<-computeChIPProfile(occupBEAFArt,chipValidationArtBEAF,cores=cores)



chipCTCF<-computeChIPProfile(occupCTCF,chipValidationCTCF,cores=cores)
chipsuhw<-computeChIPProfile(occupsuhw,chipValidationsuhw,cores=cores)

gofBEAF<-profileAccuracyEstimate(chipBEAF,chipValidationBEAF,cores=cores)
gofBEAFArt<-profileAccuracyEstimate(chipBEAFArt,chipValidationArtBEAF,cores=cores)

gofCTCF<-profileAccuracyEstimate(chipCTCF,chipValidationCTCF,cores=cores)
gofsuhw<-profileAccuracyEstimate(chipsuhw,chipValidationsuhw,cores=cores)


pdf("BEAF_validation.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(chipBEAFArt,chipValidationArtBEAF,chromatinState=AccessBEAF)
dev.off()
pdf("CTCF_validation.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(chipCTCF,chipValidationCTCF,chromatinState=AccessCTCF)
dev.off()

pdf("suhw_validation.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(chipsuhw,chipValidationsuhw,chromatinState=Accesssuhw)
dev.off()


#### plot selection

train <- c(3,8)
val<- c(2,10)
valArt<-c(8,9)

## extract data

trainChIPscores <- scores(chipTrainBEAF)[train]
trainChIPLoci <- loci(chipTrainBEAF)[train]
trainPred <- subOptiBEAF$ChIPProfiles[[1]][train]


validationChIPscore<- scores(chipValidationBEAF)[val]
validationChIPLoci<- loci(chipValidationBEAF)[val]
validationPred <- profiles(chipBEAF)[[1]][val]

validationChIPscoreArt<- scores(chipValidationArtBEAF)[valArt]
validationChIPLociArt<- loci(chipValidationArtBEAF)[valArt]
validationPredArt <- profiles(chipBEAFArt)[[1]][valArt]



metricsVal<-list("AUC"=NULL,"spearman"=NULL,"recall"=NULL,"MSE"=NULL)

for(i in seq_along(metricsVal)){
    bufferVal <- profiles(gofBEAF)[[1]]
    bufferTrain <- subOptiBEAF$goodnessOfFit[[1]]
    bufferArt <-profiles(gofBEAFArt)[[1]]
    metricsVal[[i]]<-list("Training chr3R"=sapply(bufferTrain,function(x, met){return(x[met])},names(metricsVal)[i]),
                          "Validation chr3R"=sapply(bufferVal,function(x, met){return(x[met])},names(metricsVal)[i]),
                          "Validation chr2R"=sapply(bufferArt,function(x, met){return(x[met])},names(metricsVal)[i]))

}


#metricsVal <-unlist(metricsVal, recursive=F)
#for(i in seq_along(metricsVal)){
#    if(grepl("MSE",names(metricsVal)[i])){

#         metricsVal[[i]]<-metricsVal[[i]]/max(metricsVal[[i]])
  #  }
#}


pdf("chromosome_withhold_setup.pdf",width=15,height=20)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,8,9,10), ncol=2, byrow=T),height=c(1.4,1.4,1.4,1.4,1.4,1.4,2,2))
par(family="sans",xpd=NA)


for(i in seq_along(trainChIPLoci)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(trainChIPLoci)[i],end(trainChIPLoci)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(trainChIPLoci)[i],(end(trainChIPLoci)[i])+1,by=5000), labels=seq(start(trainChIPLoci)[i],(end(trainChIPLoci)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1) {
  title(main="Training Profiles on chr3R",cex.main=2)
  text(x=start(trainChIPLoci)[i],y=1.05,label="A", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(trainChIPLoci[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(trainChIPLoci[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,trainChIPscores[[i]][seq(0,length(trainChIPscores[[i]]),by=100)],0)
predInd<-c(0,trainPred[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}

for(i in seq_along(validationChIPLociArt)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(validationChIPLociArt)[i],end(validationChIPLociArt)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(validationChIPLociArt)[i],(end(validationChIPLociArt)[i])+1,by=5000), labels=seq(start(validationChIPLociArt)[i],(end(validationChIPLociArt)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1){
  title(main="Validation Profiles on chr3R",cex.main=2)
  text(x=start(validationChIPLociArt)[i],y=1.05,label="B", cex=4)
  text(x=(end(validationChIPLociArt)[i])+2500,y=1.05,label="C", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(validationChIPLociArt[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(validationChIPLociArt[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,validationChIPscoreArt[[i]][seq(0,length(validationChIPscoreArt[[i]]),by=100)],0)
predInd<-c(0,validationPredArt[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}


for(i in seq_along(validationChIPLoci)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(validationChIPLoci)[i],end(validationChIPLoci)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(validationChIPLoci)[i],(end(validationChIPLoci)[i])+1,by=5000), labels=seq(start(validationChIPLoci)[i],(end(validationChIPLoci)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1){
  title(main="Validation Profiles on chr2R",cex.main=2)
  text(x=start(validationChIPLoci)[i],y=1.05,label="C", cex=4)
  text(x=(end(validationChIPLoci)[i])+2500,y=1.05,label="C", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(validationChIPLoci[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(validationChIPLoci[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,validationChIPscore[[i]][seq(0,length(validationChIPscore[[i]]),by=100)],0)
predInd<-c(0,validationPred[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}

counts<-list(c(1,2),c(3,4),c(5,6),c(7,8))
lims <- list(c(0.5,1),c(-0.5,0.85),c(0.5,1),c(0,0.05))
labs<-LETTERS[4:7]
for(i in seq_along(metricsVal)){
cols <-c("#4f9da6","#95c4c9","#ff5959","#ff9b9b","#facf5a","#fce29c","#009E73","#66c4ab")
#names(metricsVal)<-gsub("\\."," ", names(metricsVal))
#names(metricsVal)<-gsub("MSE","Norm. MSE",names(metricsVal))
par(mar=c(12,8,4,6))
#metricsVal[[i]]<-lapply(metricsVal[[i]])
par(xpd=NA)
text(x=0,y=max(lims[[i]])*1.1, labels=labs[i],cex=4)
boxplot(metricsVal[[i]],main=names(metricsVal)[i],col=cols[counts[[i]]],
frame=F,cex.axis=1.5,cex.main=1.8,ylim=lims[[i]],las=2,cex=2)
title(ylab="Scores", cex.lab=2, line=5)
#text(x=(-0.55),y=1.1, labels="C",cex=4)

}
dev.off()

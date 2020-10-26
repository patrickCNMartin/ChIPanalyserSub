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




## going to select regions from different chromosome
## we can base this on method data sets



input<-read.table("DataInputDHSgeometric.txt", sep=' ', comment.char='@', stringsAsFactors=F)


DNASequenceSet <- get(load(input[1,6]))
AccessBEAF<-get(load(input[7,7]))
AccessCTCF<-get(load(input[7,7]))
Accesssuhw<-get(load(input[7,7]))



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


pwmCTCF<-computePWMScore(gwCTCF,DNASequenceSet,chipValidationCTCF,access,cores=cores)
pwmsuhw<-computePWMScore(gwsuhe,DNASequenceSet,chipValidationsuhw,access,cores=cores)

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



#### plot selection

beafLoc <- c(5,16)
ctcfLoc <- c(2,3)
suhwLoc <- c(9,10)

pdf("thermo.pdf", width = 20, height = 9)
par(mfrow = c(3,2), xpd=NA)

thermoChIP <- list(chipValidationBEAF,chipValidationCTCF,chipValidationsuhw)
thermoPred <- list(chipBEAF,chipCTCF,chipsuhw)
thermoPWM <- list(occupBEAF,occupCTCF,occupsuhw)
Accesslist <- list(AccessBEAF,AccessCTCF,Accesssuhw)
loc <- list(beafLoc,ctcfLoc,suhwLoc)
tag <- c("BEAF-32","CTCF","Su(Hw)")

for(i in seq_along(thermoChIP)){
    cols<-c("#facf5a","#233142","#ff5959","#4f9da6","#56B4E9")
    
    for(k in seq_along(beafLoc)){
        trainChIPLoci<- loci(thermoChIP[[i]])[loc[[i]][k]]
        trainChIPscores <- scores(thermoChIP[[i]])[loc[[i]][k]]
        trainPred <- profiles(thermoPred[[i]])[[1]][loc[[i]][k]]
        trainPWM <- profiles(thermoPWM[[i]])[[1]][loc[[i]][k]]
        x<-seq(start(trainChIPLoci),end(trainChIPLoci),by=100)
        x<-c(x[1]-1,x,x[length(x)]+1)
        par(mar=c(6,3,4,1))
        plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=" ",xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
        axis(1,at=seq(start(trainChIPLoci),(end(trainChIPLoci))+1,by=5000), labels=seq(start(trainChIPLoci),(end(trainChIPLoci))+1,by=5000),cex.axis=1.6)
        title(main =tag[[i]], cex.main=2)
        #axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
        if(k ==1) {
            #title(main=tag[[i]],cex.main=2)
            text(x=start(trainChIPLoci),y=1.05,label=LETTERS[[i]], cex=4)
        }
        title(xlab=paste0("Genomic Position on : ",seqnames(trainChIPLoci)),cex.lab=1.6)
        
        
        ## plotting stuff
        noaccess<-.AccessExtract(trainChIPLoci,Accesslist[[i]])[[1]]
        for(j in seq_len(nrow(noaccess))){
            rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
            #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
            #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
        }
        
        chipInd<-c(0,trainChIPscores[[1]][seq(0,length(trainChIPscores[[1]]),by=100)],0)
        predInd<-c(0,trainPred[[1]]$ChIP,0)
        maxLocal <- maxPWMScore(thermoPWM[[i]])
        minLocal <- minPWMScore(thermoPWM[[i]])
        pwmInd <- trainPWM[[1]]$PWMScore
        pwmInd <- (pwmInd - minLocal)/ (maxLocal - minLocal)
        
        polygon(x,chipInd,col=cols[2],lwd=2)
        lines(x,predInd,col=cols[3],lwd=4)
        lines(start(trainPWM[[1]]),pwmInd, type ="h",col = "#56B4E9",lwd=1.5)
        
    }
}
dev.off()

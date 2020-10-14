################################################################################
############## Running CTCF on humans because people hard code shit ############
################################################################################

direc <- getwd()

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
library(rtracklayer)
library(MotifDb)


## sourcing scripts for analysis
source("DataHand.R")



#loading Data

DNASequenceSet<- getSeq(BSgenome.Hsapiens.UCSC.hg38)[1:24]
locus<-GRanges(names(DNASequenceSet), ranges=IRanges(1,width(DNASequenceSet)))
locus<-unlist(tile(locus, width=20000))

# filter black listing
black<-import("~/ChIP/human/hg38.blacklist.bed.gz")

setDiffOverlap<-function(query, subject){
    overlaps<-queryHits(findOverlaps(query,subject))
    idx<-seq_along(query)
    diff<-query[!idx %in% overlaps]
    return(diff)
}

locus<- setDiffOverlap(locus, black)


## chromosome splitting

chr11<- locus[which(as.character(seqnames(locus))=="chr11")]

chr18<- locus[which(as.character(seqnames(locus))=="chr18")]



pfmJASPAR<-query(query(query(MotifDb,"hsapiens"),"CTCF"),"jaspar2018")[[1]]

pfmHOMO<-query(query(query(MotifDb,"hsapiens"),"CTCF"),"HOCOMOCO")[[1]]

pfmFor <- "matrix"



### peak regions we are just going to use one type


## data from rafael but processed by encode
#peaks<-import("ChIP_data/ENCFF600CYD.bed.gz")
## import struggles so let's do it this way
peaks<-read.table("ChIP_data/ENCFF600CYD.bed.gz")
peaks <- GRanges(seqnames=peaks[,1],ranges=IRanges(peaks[,2], peaks[,3]))

## acces of combinded replicates --- Might be owrthwhile to only use one of the 2
Access <- import("Access/DNAse_hg38_astro_combindedreps.bed.gz")

## import signal of combinded signal over two replicates
signal <- import("ChIP_data/ENCFF424JNY.bigWig")




Access18 <- Access[which(as.character(seqnames(Access))=="chr18")]
Access11 <- Access[which(as.character(seqnames(Access))=="chr11")]



PO<-parameterOptions(noiseFilter="sigmoid")



stTrain <-Sys.time()
#ChIPProfileTraining <- processingChIP(signal,loci=chr18,reduce = length(chr18),parameterOptions=PO,cores=1)
ChIPProfileTraining <- processingChIP(signal,loci=chr18,reduce = 10,peaks=peaks,chromatinState = Access18,parameterOptions=PO,cores=1)
GP<- genomicProfiles(PFM=pfmJASPAR,PFMFormat=pfmFor, BPFrequency=DNASequenceSet,PWMThreshold=0.7)
optimal<-computeOptimal(GP, DNASequenceSet,ChIPProfileTraining,chromatinState=Access18,parameterOptions=PO,cores=1)
#pdf("~/ChIP/human/heatFull.pdf")
#plotOptimalHeatMaps(optimal)
#dev.off()

param <-optimal[[1]][[1]][["MSE"]]
subopti<-searchSites(optimal,param[1],param[2])
pdf("~/ChIP/human/profiles_MSE.pdf", width=14, height=4)
plotOccupancyProfile(subopti$ChIPProfiles,ChIPProfileTraining)
dev.off()
enTrain <-Sys.time()
train<-enTrain-stTrain
write(train, file="TrainTime.txt")


stVal <-Sys.time()
#ChIPProfileValidation <- processingChIP(signal,loci=chr11,reduce = length(chr11),parameterOptions=PO,cores=1)
ChIPProfileValidation <- processingChIP(signal,loci=chr11,reduce = 20,parameterOptions=PO,cores=1)





GPval<-genomicProfiles(PFM=pfmJASPAR,PFMFormat=pfmFor, BPFrequency=DNASequenceSet,boundMolecules=param[2],lambdaPWM=param[1])
gw <- computeGenomeWideScores(GPval, DNASequenceSet,Access11, cores=1)
pwm<-computePWMScore(gw, DNASequenceSet,ChIPProfileValidation, Access11, cores=1)
occup <-computeOccupancy(pwm)
chip<-computeChIPProfile(occup,ChIPProfileValidation,cores=1)
gof<-profileAccuracyEstimate(chip,ChIPProfileValidation,cores=1)
vali <- list("occup"=occup,"chip"=chip,"gof"=gof)
pdf("~/ChIP/human/validation_train10_va20.pdf")
par(mfrow=c(5,2),mar=c(2,2,2,2))
plotOccupancyProfile(chip,ChIPProfileValidation,chromatinState=Access11)
dev.off()

enVal <-Sys.time()
val<-enVal-stVal
write(val, file="TrainTime_train10_va20.txt")
## this data set seems to work well so let's just use this one then


## saving top regions for validation and training as bed files exported in anchor directory for now

save(ChIPProfileTraining, file="ChIPProfileTraining_chr18_train10_va20.Rda")
save(ChIPProfileValidation, file="ChIPProfileValidation_chr11_train10_va20.Rda")
save(optimal,file="Optimal_training_chr18_train10_va20.Rda")
save(vali,file="Validation_chr11_train10_va20.Rda")

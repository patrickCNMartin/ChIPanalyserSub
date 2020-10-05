#!/usr/bin/Rscript3.5.0

### Loading Libraries and Scripts

direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(rtracklayer)
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



#### Data loading
DNASequenceSet<-getSeq(BSgenome.Dmelanogaster.UCSC.dm6)[1:6]


pfms<-get(load("pfmsall.Rda"))

accessFiles<-dir(pattern="ATACAccess")
Access<-vector("list",length(accessFiles))

for(i in seq_along(Access)){
    Access[[i]]<-get(load(accessFiles[i]))
    Access[[i]]<-GRanges(Access[[i]])
}
names(Access)<-sapply(strsplit(accessFiles,".Rda"),"[[",1)


#peaks <- c("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/HOX/data/HOX_ChIP_Kc167_stable/")


#ocus<-get(load("HoX_Loci.Rda"))
#rm(loci)


setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)


################################################################################
################################################################################
################################################################################
## Common regions building when needed but it has been done
## see above as those are the saved regions
DNASequenceSet<-DNASequenceSet[1:6]
start <- vector("list", length(DNASequenceSet))
end <- vector("list", length(DNASequenceSet))
tileSize<-20000
localNames <-c()
for(i in seq_along(DNASequenceSet)){
    div <- length(DNASequenceSet[[i]]) %/% tileSize
    rem <- length(DNASequenceSet[[i]]) %% tileSize
    if(div!=0){
        start[[i]]<-seq(from=1,by=tileSize, length.out=div)
        end[[i]]<-c(seq(from=tileSize+1,by=tileSize,length.out=div-1),start[[i]][div]+rem)
    } else {
        start[[i]]<-1
        end[[i]]<-length(DNASequenceSet[[i]])
    }
    localNames <- c(localNames,rep(names(DNASequenceSet)[[i]],length(start[[i]])))
}

localSequence<- GRanges(seqnames=localNames,
    range=IRanges(start=unlist(start), end=(unlist(end)-1)))
names(localSequence)<-paste0(as.character(seqnames(localSequence)),":",start(localSequence),"..",end(localSequence))


## lets do some intersecting

#locus<-localSequence
#locusList<-vector("list",3)
#for(i in seq_along(peaks)){
    #locusList[[i]]<-locus[queryHits(findOverlaps(locus,peaks[[i]]))]
#    locusList[[i]]<-union(locus,peaks[[i]])
#    locusList[[i]]<-locusList[[i]][!duplicated(locusList[[i]])]
#}

#nah that might have been a mistake. So we be expanding

#locusList<-union(union(peaks[[1]], peaks[[2]]), peaks[[3]])

#locus<-locus[queryHits(findOverlaps(locus, locusList))]
#locus<-locus[!duplicated(locus)]

#save(locus,file="HoX_Loci.Rda")


#### anoher test for this shit because for some reason i cannot repordue this step100_ChIPTraining
peak <- import("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/HOX/data/HOX_ChIP_Kc167_stable/GSE122573_Hox_macs_q1e-10.bed.gz")

peak <- peak[peak$name !="Hth_stable"]
locus <- localSequence[queryHits(findOverlaps(localSequence, peak))]
locus <- locus[!duplicated(locus)]

## filtering some crap loci that are messing shit up
filterTag <-c("chr3R:9780001..9800000","chrX:15980001..16000000","chr3R:10600001..10620000","chr3R:10520001..10540000","chr3R:10560001..10580000")

locus<-locus[which(!names(locus)%in%filterTag)]
#### let's get the original signal
#because I just dont get what is happening here

signal<-c("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/HOX/data/HOX_ChIP_Kc167_stable/GSE122573_AbdB_stable.bedgraph.gz",
          "/home/pm16057/ChIPanalyser/ChIPanalyserFinal/HOX/data/HOX_ChIP_Kc167_stable/GSE122573_Ubx_stable.bedgraph.gz",
          "/home/pm16057/ChIPanalyser/ChIPanalyserFinal/HOX/data/HOX_ChIP_Kc167_stable/GSE122573_Dfd_stable.bedgraph.gz")

################################################################################
################################################################################
################################################################################


## let's get this party started
cores<-40
reduce<-3838
method<-c("AUC","MSE","recall","spearman")
for(i in seq_along(pfms)){
    print(names(pfms)[i])
    print(paste(i ,"pre chip"))
    PO<-parameterOptions(noiseFilter="sigmoid")
    chipProfile<-processingChIP(signal[i],locus,reduce=reduce,
                                  parameterOptions=PO,cores=cores)
    print(paste(i ,"post chip"))
    TrainScore<-scores(chipProfile)[1:10]
    ValidationScore<-scores(chipProfile)

    Trainloci<-loci(chipProfile)[1:10]
    Validationloci<-loci(chipProfile)


    .scores(chipProfile)<-TrainScore
    .loci(chipProfile)<-Trainloci
     TrainProfile<-chipProfile

     .scores(chipProfile)<-ValidationScore
     .loci(chipProfile)<-Validationloci
     ValidationProfile<-chipProfile

    GPPinitial<-genomicProfiles(PFM=pfms[[i]],PFMFormat="matrix",BPFrequency=DNASequenceSet,stepSize=100)
    for(j in seq_along(Access)){
      print(paste("Access Number =", names(Access)[j]))
      optimal <- computeOptimal(GPPinitial,DNASequenceSet,TrainProfile,Access[[j]],PO,cores=cores)
      print("post optimal")

        ### Computing Optimal Parameters ###
        filename<-paste0("../res/HoX_",names(pfms)[i],"_")

        save(optimal,file=paste0(filename,names(Access)[j],"_step100_OptimalOutputTraining.Rda"))
        save(chipProfile, file=paste0(filename,names(Access)[j],"step100_ChIPTraining.Rda"))


        ## optimal param for validation
        ## let assume you have more than one

        for(meth in method){

        param<-optimal[[1]][[1]][[meth]]

        lambda<-param[1]
        bm<-param[2]

        GPP<-genomicProfiles(PFM=pfms[[i]],PFMFormat="matrix",
        BPFrequency=DNASequenceSet,PWMThreshold=0.7,lambdaPWM=lambda,boundMolecules=bm)
        gw<-computeGenomeWideScores(GPP,DNASequenceSet,Access[[j]],cores=cores)

        pwm<-computePWMScore(gw,DNASequenceSet,ValidationProfile,Access[[j]],cores=cores)
        occup<-computeOccupancy(pwm)
        chip<-computeChIPProfile(occup,ValidationProfile,cores=cores)
        gof<-profileAccuracyEstimate(chip,ValidationProfile,cores=cores)
        optimalList<-list("Occupancy"=occup,"ChIPProfile"=chip,"Gof"=gof)
        print("post Validation")

        save(optimalList,file=paste0(filename,names(Access)[j],"_",meth,"OptimalOutputValidation.Rda"))
        save(chipProfile, file=paste0(filename,names(Access)[j],"_",meth,"ChIPValidation.Rda"))

    }
  }
}

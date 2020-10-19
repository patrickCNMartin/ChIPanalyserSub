################################################################################
########################## Method Comp scoring #################################
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
setwd("/home/patrickmartin/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)



### loading data

## ChIPanalyser
print("Start Loading Data")
ChIPTraining <-get(load("/home/patrickmartin/methodComp/ChIPanal/ChIPProfileTraining_chr18_train10_va20.Rda"))
ChIPValidation <-get(load("/home/patrickmartin/methodComp/ChIPanal/ChIPProfileValidation_chr11_train10_va20.Rda"))
optimal <-get(load("/home/patrickmartin/methodComp/ChIPanal/Optimal_training_chr18_train10_va20.Rda"))
validation <-get(load("/home/patrickmartin/methodComp/ChIPanal/Validation_chr11_train10_va20.Rda"))
print("ChIPanal Check")

## PIQ
PIQ <-import("~/methodComp/PIQ/astro_mergedBAM_CTCF-calls.all.bed")
print("PIQ Check")
##msCENTIPEDE
centi <- read.table("/home/patrickmartin/methodComp/CENTIPEDE/fimo_CTCF_astro_clean_msCentipede_binding_posterior.txt.gz",header=T)
centi <- GRanges(centi[,1],ranges=IRanges(centi[,2],centi[,3]))
print("CENTI check")
## Catchitt

catchitt <- read.table("/home/patrickmartin/methodComp/Catchitt/Predictions.tsv.gz")
print("Catchitt check ")



#### scoring

## over "training"
topTrain <-loci(ChIPTraining)
topVal <- loci(ChIPValidation)

topTrainScores<-scores(ChIPTraining)
topValScores <- scores(ChIPValidation)
print("Score Extraction Completed")

#### PIQ Training
print("start PIQ Training")
piqTrain <-PIQ[queryHits(findOverlaps(PIQ, topTrain))]

#piq<-piqTrain
top<-topTrainScores
topLocus<-topTrain
## generating templates
prediction <- vector("list",length(top))
for(i in seq_along(top)){
    local <-rep(0,length(top[[i]]))
    localLocus <- topLocus[i]
    localPrediction <- PIQ[queryHits(findOverlaps(PIQ, localLocus))]
    start <- start(localPrediction)-start(localLocus)
    if(any(start <0))start[start <0] <- 1
    end <- end(localPrediction)-start(localLocus)
    if(any(end > width(localLocus)))end[end> width(localLocus)] <- width(localLocus) 
    for(j in seq_along(start)){
        local[start[j]:end[j]]<-1
    }
    prediction[[i]]<-roll_mean(local,n=100,align="center")
    prediction[[i]]<-prediction[[i]][seq(1,length(prediction[[i]]),by=100)]
}



save(prediction,file="PIQ_profiles_Training_chr18_train10_va20.Rda")
print("PIQ Training Done")

### PIQ Validation
print("PIQ Validation")
piqVal <-PIQ[queryHits(findOverlaps(PIQ, topVal))]

piq<-piqVal
top<-topValScores
topLocus<-topVal
## generating templates
prediction <- vector("list",length(top))
for(i in seq_along(top)){
    local <-rep(0,length(top[[i]]))
    localLocus <- topLocus[i]
    localPrediction <- PIQ[queryHits(findOverlaps(PIQ, localLocus))]
    start <- start(localPrediction)-start(localLocus)
    if(any(start <0))start[start <0] <- 1
    end <- end(localPrediction)-start(localLocus)
    if(any(end > width(localLocus)))end[end> width(localLocus)] <- width(localLocus) 
    
    for(j in seq_along(start)){
        local[start[j]:end[j]]<-1
    }
    prediction[[i]]<-roll_mean(local,n=100,align="center")
    prediction[[i]]<-prediction[[i]][seq(1,length(prediction[[i]]),by=100)]
}

save(prediction,file="PIQ_profiles_Validation_chr11_train10_va20.Rda")
print("PIQ Validation Done")

### msCENTIPEDE
print("CENTIPEDE Training start ")

top<-topTrainScores
topLocus<-topTrain
## generating templates
prediction <- vector("list",length(top))
for(i in seq_along(top)){
    local <-rep(0,length(top[[i]]))
    localLocus <- topLocus[i]
    localPrediction <- centi[queryHits(findOverlaps(centi, localLocus))]
    start <- start(localPrediction)-start(localLocus)
    if(any(start <0))start[start <0] <- 1
    end <- end(localPrediction)-start(localLocus)
    if(any(end > width(localLocus)))end[end> width(localLocus)] <- width(localLocus) 
    for(j in seq_along(start)){
        local[start[j]:end[j]]<-1
    }
    prediction[[i]]<-roll_mean(local,n=100,align="center")
    prediction[[i]]<-prediction[[i]][seq(1,length(prediction[[i]]),by=100)]
}



save(prediction,file="msCENTIPEDE_profiles_Training_chr18_train10_va20.Rda")
print("CENTIPEDE training done ")


print("CENTI validation start ")

top<-topValScores
topLocus<-topVal
## generating templates
prediction <- vector("list",length(top))
for(i in seq_along(top)){
    local <-rep(0,length(top[[i]]))
    localLocus <- topLocus[i]
    localPrediction <- centi[queryHits(findOverlaps(centi, localLocus))]
    start <- start(localPrediction)-start(localLocus)
    if(any(start <0))start[start <0] <- 1
    end <- end(localPrediction)-start(localLocus)
    if(any(end > width(localLocus)))end[end> width(localLocus)] <- width(localLocus) 
    
    for(j in seq_along(start)){
        local[start[j]:end[j]]<-1
    }
    prediction[[i]]<-roll_mean(local,n=100,align="center")
    prediction[[i]]<-prediction[[i]][seq(1,length(prediction[[i]]),by=100)]
}

save(prediction,file="msCENTIPEDE_profiles_Validation_chr11_train10_va20.Rda")
print("CENTI validation done ")
### catchitt

catchittGR<-GRanges(catchitt[,1],ranges=IRanges(catchitt[,2],catchitt[,2]+100),pred=catchitt[,3])

catSub<- catchittGR[queryHits(findOverlaps(catchittGR,topVal))]

top<-topValScores
topLocus<-topVal
## generating templates

prediction<-vector("list", length(topVal))
for(i in seq_along(prediction)){
    local <-rep(0,length(top[[i]]))
    localLocus <- topLocus[i]
    print(i)
    localPrediction <- unique(catSub[queryHits(findOverlaps(catSub, localLocus))])
    start <-start(localPrediction)-start(localLocus)
    if(any(start <0))start[start <0] <- 1
    end <- end(localPrediction)-start(localLocus)
    if(any(end > width(localLocus)))end[end> width(localLocus)] <- width(localLocus) 
    
    for(j in seq_along(start)){
        local[start[j]:end[j]]<-localPrediction$pred[j]
    }
    #prediction[[i]]<-roll_mean(local,n=100,align="center")
    prediction[[i]]<-local[seq(1,length.out=199)]
    
}


save(prediction,file="Catchitt_profiles_Validation_chr11_train10_va20.Rda")
print("catchitt validation done ")

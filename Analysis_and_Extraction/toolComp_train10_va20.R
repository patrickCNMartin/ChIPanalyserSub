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
source("DataHand.R")



### loading data

## ChIPanalyser
print("Start Loading Data")
ChIPTraining <-get(load("/home/pm16057/methodComp/ChIPanal/ChIPProfileTraining_chr18_train10_va20.Rda"))
ChIPValidation <-get(load("/home/pm16057/methodComp/ChIPanal/ChIPProfileValidation_chr11_train10_va20.Rda"))
optimal <-get(load("/home/pm16057/methodComp/ChIPanal/Optimal_training_chr18_train10_va20.Rda"))
validation <-get(load("/home/pm16057/methodComp/ChIPanal/Validation_chr11_train10_va20.Rda"))
print("ChIPanal Check")

## PIQ
PIQ <-import("~/methodComp/PIQ/astro_mergedBAM_CTCF-calls.all.bed")
print("PIQ Check")
##msCENTIPEDE
centi <- read.table("/home/pm16057/methodComp/CENTIPEDE/fimo_CTCF_astro_clean_msCentipede_binding_posterior.txt.gz",header=T)
centi <- GRanges(centi[,1],ranges=IRanges(centi[,2],centi[,3]))
print("CENTI check")
## Catchitt

catchitt <- read.table("/home/pm16057/methodComp/Catchitt/Predictions.tsv.gz")
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
subs<-seq(min(sapply(top, min)),max(sapply(top, max)), length.out=21)
acc<-c()
rec<-c()
auc<-c()
mcc<-c()
prec<-c()
fscore<-c()
mse<-c()
for(j in seq_along(prediction)){
locusProfile <- top[[j]]
locusProfile<- roll_mean(locusProfile,100,by=100)

matpred<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))
matloc<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))

for(i in seq_along(subs)){

    localProfile<-rep(0,length(locusProfile))
    ## thresh extraction

    localProfile[locusProfile>=subs[i]]<-1
    if(all(localProfile==1)){
         localProfile[1]<-0
    }
    if(all(localProfile==0)){
         localProfile[1]<-1
    }
    matpred[,i]<-factor(prediction[[j]])
    matloc[,i]<-factor(localProfile)
}

## Predictions and performance
pred<-prediction(matpred,matloc)
prec<-c(prec,mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T)))
rec<-c(rec,mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T)))
fscore<-c(fscore,mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T)))
acc<-c(acc,mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T)))
mcc<-c(mcc,mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T)))
auc<-c(auc,mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T)))
mse<-c(mse,.intMSE(prediction[[j]],locusProfile))
}

PIQlist <- list("fullpred"=pred,"mse"=mse,"auc"=auc,"rec"=rec,"fscore"=fscore,"acc"=acc,"mcc"=mcc,"prec"=prec,"locus"=topLocus,"scores"=top)
save(PIQlist,file="PIQ_scores_Training_chr18_train10_va20.Rda")
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
subs<-seq(min(sapply(top, min)),max(sapply(top, max)), length.out=21)
acc<-c()
rec<-c()
auc<-c()
mcc<-c()
prec<-c()
fscore<-c()
mse<-c()
for(j in seq_along(prediction)){
locusProfile <- top[[j]]
locusProfile<- roll_mean(locusProfile,100,by=100)
matpred<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))
matloc<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))

for(i in seq_along(subs)){

    localProfile<-rep(0,length(locusProfile))
    ## thresh extraction

    localProfile[locusProfile>=subs[i]]<-1
    if(all(localProfile==1)){
         localProfile[1]<-0
    }
    if(all(localProfile==0)){
         localProfile[1]<-1
    }
    matpred[,i]<-factor(prediction[[j]])
    matloc[,i]<-factor(localProfile)
}

## Predictions and performance
pred<-prediction(matpred,matloc)
prec<-c(prec,mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T)))
rec<-c(rec,mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T)))
fscore<-c(fscore,mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T)))
acc<-c(acc,mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T)))
mcc<-c(mcc,mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T)))
auc<-c(auc,mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T)))
mse<-c(mse,.intMSE(prediction[[j]],locusProfile))
}

PIQlist <- list("fullpred"=pred,"mse"=mse,"auc"=auc,"rec"=rec,"fscore"=fscore,"acc"=acc,"mcc"=mcc,"prec"=prec,"locus"=topLocus,"scores"=top)
save(PIQlist,file="PIQ_scores_Validation_chr11_train10_va20.Rda")
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
subs<-seq(min(sapply(top, min)),max(sapply(top, max)), length.out=21)
acc<-c()
rec<-c()
auc<-c()
mcc<-c()
prec<-c()
fscore<-c()
mse<-c()
for(j in seq_along(prediction)){
locusProfile <- top[[j]]
locusProfile<- roll_mean(locusProfile,100,by=100)
matpred<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))
matloc<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))

for(i in seq_along(subs)){

    localProfile<-rep(0,length(locusProfile))
    ## thresh extraction

    localProfile[locusProfile>=subs[i]]<-1
    if(all(localProfile==1)){
         localProfile[1]<-0
    }
    if(all(localProfile==0)){
         localProfile[1]<-1
    }
    matpred[,i]<-factor(prediction[[j]])
    matloc[,i]<-factor(localProfile)
}

## Predictions and performance
pred<-prediction(matpred,matloc)
prec<-c(prec,mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T)))
rec<-c(rec,mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T)))
fscore<-c(fscore,mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T)))
acc<-c(acc,mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T)))
mcc<-c(mcc,mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T)))
auc<-c(auc,mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T)))
mse<-c(mse,.intMSE(prediction[[j]],locusProfile))
}

CENTIList <- list("fullpred"=pred,"mse"=mse,"auc"=auc,"rec"=rec,"fscore"=fscore,"acc"=acc,"mcc"=mcc,"prec"=prec,"locus"=topLocus,"scores"=top)
save(CENTIList,file="msCENTIPEDE_scores_Training_chr18_train10_va20.Rda")
print("CENTIPEDE training done ")

### PIQ Validation
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
subs<-seq(min(sapply(top, min)),max(sapply(top, max)), length.out=21)
acc<-c()
rec<-c()
auc<-c()
mcc<-c()
prec<-c()
fscore<-c()
mse<-c()
for(j in seq_along(prediction)){
locusProfile <- top[[j]]
locusProfile<- roll_mean(locusProfile,100,by=100)
matpred<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))
matloc<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))

for(i in seq_along(subs)){

    localProfile<-rep(0,length(locusProfile))
    ## thresh extraction

    localProfile[locusProfile>=subs[i]]<-1
    if(all(localProfile==1)){
         localProfile[1]<-0
    }
    if(all(localProfile==0)){
         localProfile[1]<-1
    }
    matpred[,i]<-factor(prediction[[j]])
    matloc[,i]<-factor(localProfile)
}

## Predictions and performance
pred<-prediction(matpred,matloc)
prec<-c(prec,mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T)))
rec<-c(rec,mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T)))
fscore<-c(fscore,mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T)))
acc<-c(acc,mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T)))
mcc<-c(mcc,mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T)))
auc<-c(auc,mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T)))
mse<-c(mse,.intMSE(prediction[[j]],locusProfile))
}

CENTIList <- list("fullpred"=pred,"mse"=mse,"auc"=auc,"rec"=rec,"fscore"=fscore,"acc"=acc,"mcc"=mcc,"prec"=prec,"locus"=topLocus,"scores"=top)
save(CENTIList,file="msCENTIPEDE_scores_Validation_chr11_train10_va20.Rda")
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


subs<-seq(min(sapply(top, min)),max(sapply(top, max)), length.out=21)
acc<-c()
rec<-c()
auc<-c()
mcc<-c()
prec<-c()
fscore<-c()
mse<-c()
for(j in seq_along(prediction)){
locusProfile <- top[[j]]
locusProfile<- roll_mean(locusProfile,100,align="center")
locusProfile <- locusProfile[seq(1,length(locusProfile),by=100)]
matpred<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))
matloc<-matrix(0,ncol=length(subs), nrow=length(prediction[[j]]))

for(i in seq_along(subs)){

    localProfile<-rep(0,length(locusProfile))
    localPrediction<-rep(0,length(prediction[[j]]))
    ## thresh extraction

    localProfile[locusProfile>=subs[i]]<-1
    localPrediction[prediction[[j]]>=subs[i]]<-1

    if(all(localProfile==1)){
         localProfile[1]<-0
    }
    if(all(localProfile==0)){
         localProfile[1]<-1
    }
    matpred[,i]<-factor(localPrediction)
    matloc[,i]<-factor(localProfile)
}

## Predictions and performance
pred<-prediction(matpred,matloc)
prec<-c(prec,mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T)))
rec<-c(rec,mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T)))
fscore<-c(fscore,mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T)))
acc<-c(acc,mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T)))
mcc<-c(mcc,mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T)))
auc<-c(auc,mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T)))
mse<-c(mse,.intMSE(prediction[[j]],locusProfile))
}

Catshit <- list("fullpred"=pred,"mse"=mse,"auc"=auc,"rec"=rec,"fscore"=fscore,"acc"=acc,"mcc"=mcc,"prec"=prec,"locus"=topLocus,"scores"=top)
save(Catshit,file="Catchitt_scores_Validation_chr11_train10_va20.Rda")
print("catchitt validation done ")

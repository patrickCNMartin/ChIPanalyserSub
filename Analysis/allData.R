################################################################################
###################### Lets get all this data in R !! ##########################
################################################################################

# Data set up starting with the usual

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




## custom loading function to use in apralell

.splitDNARanges <- function(DNASequenceSet,cores){
    SplitSeq <- floor(seq(1,length(DNASequenceSet),length.out=cores+1))
    splitDNASequenceSet<-vector("list",(cores))
    start <- SplitSeq[1:(length(SplitSeq)-1)]
    end <- c(SplitSeq[2:(length(SplitSeq)-1)]-1,SplitSeq[length(SplitSeq)])

    for(i in seq_len(cores)){
        splitDNASequenceSet[[i]]<-DNASequenceSet[start[i]:end[i]]
    }
    return(splitDNASequenceSet)
}

dataLoadOld<-function(directory,method=c("AUCMean","MSEMean")){
   sets<-vector("list",length(directory))
   for(i in seq_along(directory)){
      sets[[i]]<-vector("list",length(method))
      buffer<-get(load(directory[i]))
      for(j in seq_along(method)){
          sets[[i]][[j]]<-as.vector(buffer[[1]][[2]][[method[j]]])
      }
      names(sets[[i]])<-method

   }
   names(sets)<-sapply(strsplit(directory,"sigmoid"),"[[",1)
   return(sets)
}





setwd("peakWindowreduce20")


## DHS data loading for
#filesChIP<-dir()[grepl("top10opti3293", dir()) & !grepl("cont", dir()) & grepl("ChIP", dir()) & grepl("AUC", dir())]
filesOptimal<-dir()[grepl("top10opti3293", dir()) & !grepl("cont", dir()) & !grepl("ChIP", dir())& grepl("MSE", dir())]
filesTrain<-dir()[grepl("top10opti3293", dir()) & !grepl("cont", dir()) & !grepl("ChIP", dir())& grepl("Training", dir())]
aucMatDHS<-vector("list", length(filesOptimal))
aucMatDHSTrain<-vector("list", length(filesTrain))
names(aucMatDHS)<-sapply(strsplit(filesOptimal,"reduce"),"[[",1)
names(aucMatDHSTrain)<-sapply(strsplit(filesOptimal,"reduce"),"[[",1)

for(i in seq_along(filesOptimal)){
    buf<-get(load(filesOptimal[i]))
    buf<-profiles(buf$Gof)[[1]]
    aucMatDHS[[i]]<-sapply(buf,"[[","spearman")

  #  buf<-get(load(filesTrain[i]))
  #  buf<-as.vector(buf[[1]][[2]][["MSE"]])
  #  aucMatDHSTrain[[i]]<-buf
}


for(i in seq_along(filesTrain)){


    buf<-get(load(filesTrain[i]))
    buf<-as.vector(buf[[1]][[2]][["spearman"]])
    aucMatDHSTrain[[i]]<-buf
}


aucMatDHS<-do.call("cbind",aucMatDHS)
aucMatDHSTrain<-do.call("cbind",aucMatDHSTrain)


filesOptimal<-dir()[grepl("top10opti3293", dir()) & grepl("cont", dir()) & !grepl("ChIP", dir())& grepl("MSE", dir())]
filesTrain<-dir()[grepl("top10opti3293", dir()) & grepl("cont", dir()) & !grepl("ChIP", dir())& grepl("Training", dir())]

aucMatcont<-vector("list", length(filesOptimal))
names(aucMatcont)<-sapply(strsplit(filesOptimal,"reduce"),"[[",1)
aucMatcontTrain<-vector("list", length(filesTrain))
names(aucMatcontTrain)<-sapply(strsplit(filesTrain,"reduce"),"[[",1)
for(i in seq_along(filesOptimal)){
    buf<-get(load(filesOptimal[i]))
    buf<-profiles(buf$Gof)[[1]]
    aucMatcont[[i]]<-sapply(buf,"[[","spearman")
    #buf<-get(load(filesTrain[i]))
  #  buf<-as.vector(buf[[1]][[2]][["MSE"]])
    #aucMatcontTrain[[i]]<-buf
}

for(i in seq_along(filesTrain)){
    buf<-get(load(filesTrain[i]))
    buf<-as.vector(buf[[1]][[2]][["spearman"]])
    aucMatcontTrain[[i]]<-buf
}


aucMatcont<-do.call("cbind",aucMatcont)
aucMatcontTrain<-do.call("cbind",aucMatcontTrain)


setwd("../peakNULLreduce20")

filesOptimal<-dir()[grepl("top10opti3293", dir()) & grepl("NULL", dir()) & !grepl("ChIP", dir())& grepl("MSE", dir())]
filesTrain<-dir()[grepl("top10opti3293", dir()) & grepl("NULL", dir()) & !grepl("ChIP", dir())& grepl("Training", dir())]

aucMatNULL<-vector("list", length(filesOptimal))
names(aucMatNULL)<-sapply(strsplit(filesOptimal,"reduce"),"[[",1)
aucMatNULLTrain<-vector("list", length(filesTrain))
names(aucMatNULLTrain)<-sapply(strsplit(filesTrain,"reduce"),"[[",1)
for(i in seq_along(filesOptimal)){
    buf<-get(load(filesOptimal[i]))
    buf<-profiles(buf$Gof)[[1]]
    aucMatNULL[[i]]<-sapply(buf,"[[","spearman")
  #  buf<-get(load(filesTrain[i]))
    #buf<-as.vector(buf[[1]][[2]][["MSE"]])
  #  aucMatNULLTrain[[i]]<-buf
}


for(i in seq_along(filesTrain)){

    buf<-get(load(filesTrain[i]))
    buf<-as.vector(buf[[1]][[2]][["spearman"]])
    aucMatNULLTrain[[i]]<-buf
}


aucMatNULL<-do.call("cbind",aucMatNULL)
aucMatNULLTrain<-do.call("cbind",aucMatNULLTrain)

setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/")
save(aucMatDHS, file="spearmanDHS_validated.Rda")
save(aucMatcont, file="spearmancont_validated.Rda")
save(aucMatNULL, file="spearmanNULL_validated.Rda")

setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/")
save(aucMatDHSTrain, file="spearmanDHS_training.Rda")
save(aucMatcontTrain, file="spearmancont_training.Rda")
save(aucMatNULLTrain, file="spearmanNULL_training.Rda")








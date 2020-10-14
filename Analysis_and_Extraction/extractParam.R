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

extractParamOld<-function(optimal,allregions=TRUE,
                       bm=c(1,10,20,50,100,200,500, 1000, 2000,5000,10000,20000, 50000, 100000,200000,500000, 1000000),
                       lambda=seq(0.25,5, by=0.25),
                       method="MSE"){



    if(allregions){
       buffername<-sapply(strsplit(colnames(optimal[[1]]),"reduce"),"[[",1)

       buffername<-gsub("[\\(\\)]", "", buffername)
       colnames(optimal[[1]])<-gsub("[\\(\\)]", "",  colnames(optimal[[1]]))
       colnames(optimal[[2]])<-gsub("[\\(\\)]", "",  colnames(optimal[[2]]))

       buffername<-unique(buffername)



       optimalParam<-vector("list",length(buffername))
       for(i in seq_along(buffername)){

           dhs<-optimal[[1]][,grep(buffername[i],colnames(optimal[[1]]))]
           null<-optimal[[2]][,grep(buffername[i],colnames(optimal[[2]]))]

           dhstop<-apply(dhs,2, function(x,method){
                             if(method %in% c("KsDist","MSE","geometric")){
                                 dim(x)<-c(length(lambda),length(bm))
                                 local<-which(x==min(x),arr.ind=TRUE)
                                 if(nrow(local)>1){
                                   local<-local[nrow(local),]
                                 }
                                 score<-x[local[1],local[2]]

                             } else{
                               dim(x)<-c(length(lambda),length(bm))
                               local<-which(x==max(x),arr.ind=TRUE)
                               if(nrow(local)>1){
                                 local<-local[nrow(local),]
                               }
                              score<-x[local[1],local[2]]
                             }
                            return(list("idx"=local,"score"=score))},method=method)

          ScoreOnly <-sapply(dhstop, function(x){
                                 x$score})



          if(method %in% c("KsDist","MSE","geometric")){
          region<-names(dhstop)[which(ScoreOnly==min(ScoreOnly))]
          region<-strsplit(region,"reduce")[[1]][2]
          dhstop<-dhstop[[which(ScoreOnly==min(ScoreOnly))]]

        }else{
          region<-names(dhstop)[which(ScoreOnly==max(ScoreOnly))]
          region<-strsplit(region,"reduce")[[1]][2]
          dhstop<-dhstop[[which(ScoreOnly==max(ScoreOnly))]]
        }

          boundMolecules<-bm[dhstop[[1]][2]]
          lambdas<-lambda[dhstop[[1]][1]]
          score<-dhstop[[2]]

          bufferOptimal<-list("bm"=boundMolecules,"lambda"=lambdas,"score"=score,"region"=region)

          nulltop<-apply(null,2, function(x,method){
                            if(method %in% c("KsDist","MSE","geometric")){
                                dim(x)<-c(length(lambda),length(bm))
                                local<-which(x==min(x),arr.ind=TRUE)
                                if(nrow(local)>1){
                                  local<-local[nrow(local),]
                                }
                                  score<-x[local[1],local[2]]
                            } else{
                              dim(x)<-c(length(lambda),length(bm))
                              local<-which(x==max(x),arr.ind=TRUE)
                              if(nrow(local)>1){
                                local<-local[nrow(local),]
                              }
                                score<-x[local[1],local[2]]
                            }
                           return(list("idx"=local,"score"=score))},method=method)
       ScoreOnly <-sapply(nulltop, function(x){
                                    x$score})
       if(method %in% c("KsDist","MSE","geometric")){
         region<-names(nulltop)[which(ScoreOnly==min(ScoreOnly))]
         region<-strsplit(region,"reduce")[[1]][2]
         nulltop<-nulltop[[which(ScoreOnly==min(ScoreOnly))]]

       }else{
         region<-names(nulltop)[which(ScoreOnly==min(ScoreOnly))]
         region<-strsplit(region,"reduce")[[1]][2]
         nulltop<-nulltop[[which(ScoreOnly==max(ScoreOnly))]]
       }


         boundMolecules<-bm[nulltop[[1]][2]]
         lambdas<-lambda[nulltop[[1]][1]]
         score<-nulltop[[2]]
         bufferOptimalnull<-list("bm"=boundMolecules,"lambda"=lambdas,"score"=score,"region"=region)
         optimalParam[[i]]<-list("dhs"=bufferOptimal,"null"=bufferOptimalnull)


       }
       names(optimalParam)<-buffername


    } else{
      buffername<-colnames(optimal[[1]])


      optimalParam<-vector("list",length(buffername))
      names(optimalParam)<-buffername
      for(i in seq_along(buffername)){
          dhs<-optimal[[1]][,i]
          null<-optimal[[2]][,i]

        if(method %in% c("KsDist","MSE","geometric")){
                    dim(dhs)<-c(length(lambda),length(bm))
                    localdhs<-which(dhs==min(dhs),arr.ind=TRUE)
                    if(nrow(localdhs)>1){
                      localdhs<-localdhs[nrow(localdhs),]
                    }
                      scoredhs<-x[localdhs[1],localdhs[2]]
                    dim(null)<-c(length(lambda),length(bm))
                    localnull<-which(null==min(null),arr.ind=TRUE)
                    if(nrow(localnull)>1){
                      localnull<-localnull[nrow(localnull),]
                    }
                    scorenull<-null[localnull[1],localnull[2]]
        } else{
          print(i)
                  dim(dhs)<-c(length(lambda),length(bm))
                  localdhs<-which(dhs==max(dhs),arr.ind=TRUE)
                  if(nrow(localdhs)>1){
                    localdhs<-localdhs[nrow(localdhs),]
                  }
                  scoredhs<-dhs[localdhs[1],localdhs[2]]

                  dim(null)<-c(length(lambda),length(bm))
                  localnull<-which(null==max(null),arr.ind=TRUE)
                  if(nrow(localnull)>1){
                    localnull<-localnull[nrow(localnull),]
                  }
                  scorenull<-null[localnull[1],localnull[2]]
        }
         boundMolecules<-bm[localdhs[2]]
         lambdas<-lambda[localdhs[1]]
         score<-scoredhs

         bufferOptimal<-list("bm"=boundMolecules,"lambda"=lambdas,"score"=score)

         boundMolecules<-bm[localnull[2]]
         lambdas<-lambda[localnull[1]]
         score<-scorenull

         bufferOptimalnull<-list("bm"=boundMolecules,"lambda"=lambdas,"score"=score)
         optimalParam[[i]]<-list("dhs"=bufferOptimal,"null"=bufferOptimalnull)
    }
  }
    return(optimalParam)
}





################################################################################
################################################################################


extractParam<-function(data,
                      bm=c(1,10,20,50,100,200,500, 1000, 2000,5000,10000,20000, 50000, 100000,200000,500000, 1000000),
                      lambda=seq(0.25,5, by=0.25),
                      method="MSE"){
      #buildTemplate

      temp <- vector("list", ncol(data))
      names(temp)<-colnames(data)
      for(i in seq_along(temp)){
         buffer<-as.vector(data[,i])
         dim(buffer)<-c(length(lambda), length(bm))
         colnames(buffer)<-bm
         rownames(buffer)<-lambda
         #temp[[i]]<-buffer
         if(method %in% c("MSE","ks","geometric")){
             best <- which(buffer==min(buffer), arr.ind=TRUE)
             if(nrow(best)>1)best<-matrix(best[1,],ncol=2)
             temp[[i]]<-c(buffer[best],bm[best[,2]], lambda[best[,1]])
             names(temp[[i]])<-c("Score","bm","lambda")
         } else{
             best <- which(buffer==max(buffer), arr.ind=TRUE)
             if(nrow(best)>1)best<-matrix(best[1,],ncol=2)
             temp[[i]]<-c(buffer[best],bm[best[,2]], lambda[best[,1]])
             names(temp[[i]])<-c("Score","bm","lambda")
         }
      }
    return(temp)

}

latexFormat <- function(param, filename){
    ## we will expect to collapse all type of access
    filename<-paste0(filename,".txt")
    files<-gsub("_"," ",names(param[[1]]))
    #files<-gsub("modEncode","",files)
    cat("",file=filename)
    for(i in seq_along(param[[1]])){
      both<-paste0(files[i],"&",param[[1]][[i]]["bm"],"&",param[[1]][[i]]["lambda"],"&",
      round(param[[1]][[i]]["Score"],digits=3),"&",
      param[[2]][[i]]["bm"],"&",param[[2]][[i]]["lambda"],"&",
      round(param[[2]][[i]]["Score"],digits=3),"&",
      param[[3]][[i]]["bm"],"&",param[[3]][[i]]["lambda"],"&",
      round(param[[3]][[i]]["Score"],digits=3),"\\","\\","\n","\\hline","\n")

      cat(both,file=filename,append=TRUE)
    }
}


mseMatDHSTrain<-get(load("MSEDHS_training.Rda"))
mseMatcontTrain<-get(load("MSEcont_training.Rda"))
mseMatNULLTrain<-get(load("MSENULL_training.Rda"))

mse<-list(extractParam(ordered(mseMatNULLTrain)),extractParam(ordered(mseMatcontTrain)),extractParam(ordered(mseMatDHSTrain)))
latexFormat(mse,"MSE_table")


aucMatDHSTrain<-get(load("AUCDHS_training.Rda"))
aucMatcontTrain<-get(load("AUCcont_training.Rda"))
aucMatNULLTrain<-get(load("AUCNULL_training.Rda"))
auc<-list(extractParam(ordered(aucMatNULLTrain),method="AUC"),extractParam(ordered(aucMatcontTrain),method="AUC"),extractParam(ordered(aucMatDHSTrain),method="AUC"))
latexFormat(auc,"AUC_table")

recallMatDHSTrain<-get(load("recallDHS_training.Rda"))
recallMatcontTrain<-get(load("recallcont_training.Rda"))
recallMatNULLTrain<-get(load("recallNULL_training.Rda"))
recall<-list(extractParam(ordered(recallMatNULLTrain),method="recall"),extractParam(ordered(recallMatcontTrain),method="recall"),extractParam(ordered(recallMatDHSTrain),method="recall"))
latexFormat(recall,"recall_table")

spearMatDHSTrain<-get(load("spearDHS_training.Rda"))
spearMatcontTrain<-get(load("spearcont_training.Rda"))
spearMatNULLTrain<-get(load("spearNULL_training.Rda"))
spear<-list(extractParam(ordered(spearMatNULLTrain),method="spearman"),extractParam(ordered(spearMatcontTrain),method="spearman"),extractParam(ordered(spearMatDHSTrain),method="spearman"))
latexFormat(spear,"spear_table")

################################################################################
######################## Combining dm6 data sets together  #####################
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





## loading input

input<-read.table("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputDHSgeometric.txt",
                  sep=' ', comment.char='@', stringsAsFactors=F)



## Functions

splitThis<-function(data, tf=c("CTCF","Hw","BEAF"),cell=c("BG3","Kc167","S2")){
    tfs<-vector("list", length(tf))
    names(tfs)<-tf

    for(i in seq_along(tf)){
      tfs[[i]]<-vector("list",length(cell))
      names(tfs[[i]])<-cell
      local<-data[grep(tf[i],data[,2]),]
      for(j in seq_along(cell)){
          tfs[[i]][[j]]<-as.vector(as.matrix(local[grep(cell[j],local[,1]),5]))
      }
    }
    return(tfs)
}

addingChips<-function(chip,locus,cores=1){

     for(i in seq_along(chips)){
       print(names(chips)[i])
       for(j in seq_along(chips[[i]])){
          print(names(chips[[i]])[j])
          local<-vector("list", length(chips[[i]][[j]]))
         for(k in seq_along(chips[[i]][[j]])){
              print(paste("Data set",k))
              local[[k]]<-processingChIP(chips[[i]][[j]][[k]], locus, reduce=3293,cores=cores)
         }


         localScores<-scores(local[[1]])
         for(m in 2:length(local)){
             idx<-seq_along(scores(local[[m]]))
             localScores<-mapply(function(sc1,sc2,idx){
               if(length(sc1)!=length(sc2)){print(paste("length Mismatch at ",idx))
                      print(paste("length Diff",length(sc1)-length(sc2)))
                }
               return(sc1+sc2)},localScores,scores(local[[m]]),idx)
         }

         print("Added scores")
         localMax<-max(sapply(localScores,max))
         #localMin<-min(sapply(localScores,min))
         localMean<-mean(sapply(localScores, mean))

         print("metrics averaged")
            localScores<-lapply(localScores,"/",localMax)
            maxSc<-max(sapply(localScores,max))
            backGr<-mean(sapply(localScores, mean))
          buffer<-local[[1]]
          .scores(buffer)<-localScores
          maxSignal(buffer)<-maxSc
          backgroundSignal(buffer)<-backGr
          print("ChIPscore rebuild")
          chips[[i]][[j]]<-buffer

       }
     }

  return(chips)
}




## doing some stuff
chips<- splitThisShit(splitThis)

locus<-get(load(input[1,8]))

combi<-addingChips(chips,locus,cores=20)

pfms<-unique(as.vector(as.matrix(input[,3])))

DNASequenceSet<-get(load(input[1,6]))

AccessDat<-unique(as.vector(as.matrix(input[,7])))
Access<-vector("list", length(AccessDat))
names(Access)<-sapply(strsplit(AccessDat,"/cellAccess/"),"[[",2)
for(i in seq_along(Access)){
    Access[[i]]<-get(load(AccessDat[[i]]))
}

PO<-parameterOptions(noiseFilter="sigmoid",stepSize=100)

noiseFilter<-"sigmoid"
method<-c("MSE","AUC","spearman","recall")
## Data be load let's run an test

for( i in seq_along(combi)){
   GPP<-genomicProfiles(PFM=pfms[grep(names(combi)[i],pfms)],PFMFormat="JASPAR",BPFrequency=DNASequenceSet,stepSize=100)
   for(j in seq_along(combi[[i]])){
     chipProfile<-combi[[i]][[j]]

      bufferOrd <- order(sapply(scores(chipProfile),max),decreasing=T)
     TrainScore<-scores(chipProfile)[bufferOrd]
     TrainScore<-TrainScore[1:10]

     ValidationScore<-scores(chipProfile)[bufferOrd]

     Trainloci<-loci(chipProfile)[bufferOrd]
     Trainloci<-Trainloci[1:10]
     Validationloci<-loci(chipProfile)[bufferOrd]


     .scores(chipProfile)<-TrainScore
     .loci(chipProfile)<-Trainloci
     optimal <- computeOptimal(GPP,DNASequenceSet,chipProfile,Access[[grep(names(combi[[i]])[j],names(Access))]],PO,cores=cores)


       ### Computing Optimal Parameters ###
      filename <-paste0(names(combi)[i],"_", names(combi[[i]])[j],"_combinded_")
      print(filename)
       if(!is.null(filename)){
           save(optimal,file=paste0(filename,noiseFilter,"_OptimalOutputTraining.Rda"))
           save(chipProfile, file=paste0(filename,noiseFilter,"_ChIPTraining.Rda"))
       } else {
           save(optimal,file=paste0(noiseFilter,"optimalOutputTraining.Rda"))
           save(chipProfile, file=paste0(noiseFilter,"ChIPTraining.Rda"))
       }

       ## optimal param for validation
       ## let assume you have more than one

       for(meth in method){

       param<-optimal[[1]][[1]][[meth]]

       lambda<-param[1]
       bm<-param[2]

       GPP<-genomicProfiles(PFM=pfms[grep(names(combi)[i],pfms)],PFMFormat="JASPAR",
       BPFrequency=DNASequenceSet,stepSize=100,lambdaPWM=lambda, boundMolecules=bm)
       gw<-computeGenomeWideScores(GPP,DNASequenceSet,Access[[grep(names(combi[[i]])[j],names(Access))]],cores=cores)
       .scores(chipProfile)<-ValidationScore
       .loci(chipProfile)<-Validationloci
       pwm<-computePWMScore(gw,DNASequenceSet,chipProfile,Access[[grep(names(combi[[i]])[j],names(Access))]],cores=cores)
       occup<-computeOccupancy(pwm)
       chip<-computeChIPProfile(occup,chipProfile,cores=cores)
       gof<-profileAccuracyEstimate(chip,chipProfile,cores=cores)
       optimalList<-list("Occupancy"=occup,"ChIPProfile"=chip,"Gof"=gof)

       if(!is.null(filename)){
           save(optimalList,file=paste0(filename,noiseFilter,"_",meth,"OptimalOutputValidation.Rda"))
           save(chipProfile, file=paste0(filename,noiseFilter,"_",meth,"ChIPValidation.Rda"))
       } else {
           save(optimalList,file=paste0(noiseFilter,"optimalOutputValidation.Rda"))
           save(chipProfile, file=paste0(noiseFilter,"ChIPValidation.Rda"))
       }

     }

   }

}

###############################################################################
######################### QDA DNase ###########################################
###############################################################################


direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)


density <-dir()[grep("density.wig",dir())]
filename<-sapply(strsplit(density,"\\."),"[[",1)
chain<-import("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/poster_+_paper_+_thesis/objects/dm3ToDm6.over.chain")
#quants <- c(seq(0,0.9,by=0.1),0.95,0.99,0.999,0.85)
quants <- seq(0.9,0.99,l=15)
for(j in seq_along(density)){
   buffer<-read.table(density[j],skip=1,stringsAsFactors=F)
   buffer<-GRanges(seqnames=paste0("chr",buffer[,1]),
                   ranges=IRanges(buffer[,2],buffer[,3]),
                   scores=buffer[,4])
    for(i in quants){

        local<-reduce(buffer[buffer$scores>quantile(buffer$scores,i)])
        local<-unlist(liftOver(local,chain))
        save(local,file=paste0(filename[j],"_",i,".Rda"))
    }

}


## normalised scores without log

for(j in seq_along(density)){
  buffer<-read.table(density[j],skip=1,stringsAsFactors=F)
  buffer[,4]<-(buffer[,4]-min(buffer[,4]))/(max(buffer[,4])-min(buffer[,4]))
  buffer<-GRanges(seqnames=paste0("chr",buffer[,1]),
                  ranges=IRanges(buffer[,2],buffer[,3]),
                  scores=buffer[,4])
  buffer<-unlist(liftOver(buffer,chain))
  save(buffer,file=paste0(filename[j],"_MinMax_Norm_density_Scores.Rda"))
}

### max div without log


for(j in seq_along(density)){
  buffer<-read.table(density[j],skip=1,stringsAsFactors=F)
  buffer[,4]<-(buffer[,4]/max(buffer[,4]))
  buffer<-GRanges(seqnames=paste0("chr",buffer[,1]),
                  ranges=IRanges(buffer[,2],buffer[,3]),
                  scores=buffer[,4])
  buffer<-unlist(liftOver(buffer,chain))
  save(buffer,file=paste0(filename[j],"_Max_Norm_density_Scores.Rda"))
}

### norm with log

for(j in seq_along(density)){
  buffer<-read.table(density[j],skip=1,stringsAsFactors=F)

  buffer[,4]<-(buffer[,4]-min(buffer[,4]))/(max(buffer[,4])-min(buffer[,4]))
  buffer[,4]<-log2(buffer[,4]+1)
  buffer<-GRanges(seqnames=paste0("chr",buffer[,1]),
                  ranges=IRanges(buffer[,2],buffer[,3]),
                  scores=buffer[,4])
  buffer<-unlist(liftOver(buffer,chain))
  save(buffer,file=paste0(filename[j],"_Log2_MinMax_Norm_density_Scores.Rda"))
}

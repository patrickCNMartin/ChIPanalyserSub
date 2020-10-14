###############################################################################
######################### Peak overlap ########################################
###############################################################################

# loading stuff

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)

library(ChIPanalyser)

source("DataHand.R")


## Data


###################################################################################################################
###################################################################################################################
###################################################################################################################
# some Functions

addOrderIdx<-function(set,intervals){
    intervals<-rev(intervals)
    if(length(set)<intervals[1]){
        intervals<-c(length(set),intervals[2:length(intervals)])
    } else if(length(set)==intervals[1]){
         intervals<-intervals
    } else {
         intervals<-c(length(set),intervals)
    }
    ##first set
    set$order<-rep(0, length(set))

    for(i in seq_along(intervals)){
        mcols(set)[1:intervals[i],"order"]<-intervals[i]
    }
    return(set)
}


# custom function time
setDiffOverlap<-function(query, subject){
  overlaps<-queryHits(findOverlaps(query,subject))
  idx<-seq_along(query)
  diff<-query[!idx %in% overlaps]
  return(diff)
}

###################################################################################################################
###################################################################################################################
###################################################################################################################
### yet another method BUT THE FINAL ONE DO THIS ONE IF ANY
##

input<-read.table("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputDHSgeometric.txt",
                  sep=' ', comment.char='@', stringsAsFactors=F)


peaks<-input[which(!duplicated(input[,10])),]
peaks<-split(peaks[,c(1:2,10)],peaks$V1)
peaks<-lapply(peaks,function(x,tfs){
  buf<-vector("list", length(tfs))
  names(buf)<-tfs
  
  for(i in seq_along(tfs)){
    buf[[i]]<-x[grepl(tfs[i],x[,2]),3]
  }
  return(buf)
},tfs=c("CTCF","BEAF","Hw"))

DNASequenceSet<-getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
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


full<-vector("list", 3)
for(i in seq_along(peaks)){
    full[[i]]<-peaks[[i]][[1]]
    for(j in 2:3){
        full[[i]]<-union(full[[i]],peaks[[i]][[j]])
    }
    full[[i]]<-localSequence[unique(queryHits(findOverlaps(localSequence,full[[i]])))]
}

bg3<-read.table("../DNAaccess/cellAccess/BG3_DNASE_harv/BG3_Dnase_trim_005_broad_peaks.xls",header=T)
s2<-read.table("../DNAaccess/cellAccess/S2_DNASE_harv/S2_Dnase_trim_005_broad_peaks.xls",header=T)
Kc<-read.table("../DNAaccess/cellAccess/Kc_DNASE_harv/Kc_Dnase_trim_005_broad_peaks.xls",header=T)
access<-list(bg3,Kc,s2)

## checking % on accesible DNA accessibility
lapply(access, function(x){sum(x[,3]-x[,2])/1.3e8})

## checking for overlaps with accessible DNA
access<-lapply(access,function(x){GRanges(seqnames=x[,1],ranges=IRanges(x[,2],x[,3]))})
signal<-list(signal[[1]][[1]][6],signal[[1]][[2]][4], signal[[1]][[3]][5],signal[[2]][[1]][1],signal[[2]][[2]][2],
             signal[[2]][[2]][2],signal[[3]][[1]][1],signal[[3]][[2]][1],signal[[3]][[3]][3])

## ordered per cell line just in case
processed<-vector("list",9)
count<-1
for(i in seq_along(processed)){

    full[[count]]<-full[[count]][unique(queryHits(findOverlaps(full[[count]],access[[count]])))]
    processed[[i]]<-processingChIP(profile=signal[[i]],loci= full[[count]], reduce=length(full[[count]]),cores=5)
    print(i)
    if(i==3)count<-count+1 ; if(i==6)count<-count+1
    print(count )
}

## final union between different cell lines

allcell<-full[[1]]
for(i in 2:3){
    allcell<-allcell[queryHits(findOverlaps(allcell,full[[i]]))]
}
# order with same data set of loci
processed<-vector("list",9)

for(i in seq_along(signal)){
   processed[[i]]<-processingChIP(profile=signal[[i]],loci=allcell, reduce=length(allcell),cores=5)
}

intervals<-c(20,50,100,150,200,250,500,1000,2000,3000)
reorder<-vector("list", length(processed))
for(i in seq_along(processed)){
    reorder[[i]]<-addOrderIdx(loci(processed[[i]]),intervals)

}
save(reorder, file="/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/Loci_order_tagged.Rda")


## finding overlapps between top regions
top<-vector("list", length(intervals))
for(i in seq_along(top)){
    top[[i]]<-vector("list",2)

    for(j in 1:length(reorder)){
        buffer<-reorder[[j]][reorder[[j]]$order==intervals[i]]
        if(j ==1){
            top[[i]][[1]]<-buffer
            top[[i]][[2]]<-buffer
        } else{
            top[[i]][[1]]<-union(top[[i]][[1]],buffer)
            top[[i]][[2]]<-top[[i]][[2]][queryHits(findOverlaps(top[[i]][[2]],buffer))]
            print(length(top[[i]][[2]]))
        }

    }
}
###################################################################################################################
###################################################################################################################
###################################################################################################################















### another way of doing it  way of doing it
input<-read.table("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DataInputDHSgeometric.txt",
                  sep=' ', comment.char='@', stringsAsFactors=F)

#### we dont want to do this actually. we just want it for all cells the same regions
## peaks by cell line and
peaks<-input[which(!duplicated(input[,10])),]
peaks<-split(peaks[,c(1:2,10)],peaks$V1)
peaks<-lapply(peaks,function(x,tfs){
         buf<-vector("list", length(tfs))
         names(buf)<-tfs

         for(i in seq_along(tfs)){
            buf[[i]]<-x[grepl(tfs[i],x[,2]),3]
         }
         return(buf)
    },tfs=c("CTCF","BEAF","Hw"))


loci<-vector("list",9)

#processing ChIP per cell line

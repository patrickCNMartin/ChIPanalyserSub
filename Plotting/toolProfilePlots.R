direc <- getwd()

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
library(rtracklayer)
library(MotifDb)
#library(ChIPanalyser)


## sourcing scripts for analysis
setwd("/home/patrickmartin/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)


#source("DataHand.R")
## load data

catchitt <- dir()[grepl("profiles_Validation", dir()) & grepl("Catchitt", dir()) & grepl("train10", dir())]
centi <- dir()[grepl("profiles_Validation", dir()) & grepl("msCENTIPEDE", dir()) & grepl("train10", dir())]
PIQ <- dir()[grepl("profiles_Validation", dir()) & grepl("PIQ", dir()) & grepl("train10", dir())]


ChIP <- paste0("ChIPanal/","ChIPProfileValidation_chr11_train10_va20.Rda")

ChIPanal <- paste0("ChIPanal/","Validation_chr11_train10_va20.Rda")

catchitt <- get(load(catchitt))
centi <- get(load(centi))
PIQ <- get(load(PIQ))
ChIPanal <- get(load(ChIPanal))
ChIP <- get(load(ChIP))
locus <- loci(ChIP)
Access <- import("/home/patrickmartin/ChIP/human/Access/DNAse_hg38_astro_combindedreps.bed.gz")


AccessExtract<-function(subject,query){
    setLocal<-vector("list",length(subject))
    
    for(i in seq_along(subject)){
        localIntersect<-setdiff(subject[i], query)
        setLocal[[i]]<-data.frame("chr"=as.character(seqnames(localIntersect)),"start"=start(localIntersect), "end"=end(localIntersect))
    }
    names(setLocal)<-names(subject)
    return(setLocal)
}


## Extract ChIP scores
ChIPsc <- scores(ChIP)

## by name 
#loc <- "chr11:61893871..61913870"
#loci <- loci[[which(names(loci)==loc)]]
## by index 
loc <- 10
locus <- locus[loc]



## Access extract 
noaccess <- AccessExtract(locus,Access)[[1]]

cols<-c("#ff5959","#233142","#facf5a","#facf5a")
pdf("toolCompProfile.pdf",width=12,height = 14)
par(mfrow = c(4,1), mar = c(6,2,5,2),xpd=NA)


x<-seq(start(locus),end(locus),length.out=199)
x<-c(x[1]-1,x,x[length(x)]+1)
local<-ChIPsc[[loc]]
chipInd<-c(0,local[round(seq(1,length(local),length.out =199))],0)

plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(locus),end(locus)),ylim=c(0,1))

title(xlab=paste0("Genomic Position on ",as.character(seqnames(locus))),cex.lab=1.5)
title(main = "ChIPanalyser", cex.main = 3)
text(x= start(locus)-500,y=1.1,"A",cex=4)
axis(1,at=round(seq(start(locus),end(locus),length.out=10)),labels=round(seq(start(locus),end(locus),length.out=10)),cex.axis=1.5)



for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
   
}


predInd<-profiles(ChIPanal[["chip"]])[[1]][[loc]]$ChIP
predInd <- c(0,predInd[seq(1,length(predInd),length.out = 199)],0)


polygon(x,chipInd,density=NA,col="#233142",lwd=2)
lines(x,predInd,col="#ff5959",lwd=2.5)
message("chip")
#piq
plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(locus),end(locus)),ylim=c(0,1))

title(xlab=paste0("Genomic Position on ",as.character(seqnames(locus))),cex.lab=1.5)
title(main = "PIQ", cex.main = 3)
text(x= start(locus)-500,y=1.1,"B",cex=4)

axis(1,at=round(seq(start(locus),end(locus),length.out=10)),labels=round(seq(start(locus),end(locus),length.out=10)),cex.axis=1.5)



for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
    
}


predInd<-c(0,PIQ[[loc]],0)

x<-seq(start(locus),end(locus),length.out=199)
x<-c(x[1]-1,x,x[length(x)]+1)

polygon(x,chipInd,density=NA,col="#233142",lwd=2)
lines(x,predInd,col="#ff5959",lwd=2.5)
message("piq")

##CXenti 
plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(locus),end(locus)),ylim=c(0,1))

title(xlab=paste0("Genomic Position on ",as.character(seqnames(locus))),cex.lab=1.5)
title(main = "msCENTIPEDE", cex.main = 3)
text(x= start(locus)-500,y=1.1,"C",cex=4)

axis(1,at=round(seq(start(locus),end(locus),length.out=10)),labels=round(seq(start(locus),end(locus),length.out=10)),cex.axis=1.5)



for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
    
}


predInd<-c(0,centi[[loc]],0)
x<-seq(start(locus),end(locus),length.out=199)
x<-c(x[1]-1,x,x[length(x)]+1)
polygon(x,chipInd,density=NA,col="#233142",lwd=2)
lines(x,predInd,col="#ff5959",lwd=2.5)

message("centi")
### catchitt

##CXenti 
plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(locus),end(locus)),ylim=c(0,1))

title(xlab=paste0("Genomic Position on ",as.character(seqnames(locus))),cex.lab=1.5)
title(main = "Catchitt", cex.main = 3)
text(x= start(locus)-500,y=1.1,"D",cex=4)

axis(1,at=round(seq(start(locus),end(locus),length.out=10)),labels=round(seq(start(locus),end(locus),length.out=10)),cex.axis=1.5)



for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
    
}


predInd<-c(0,catchitt[[loc]],0)
x<-seq(start(locus),end(locus),length.out=199)
x<-c(x[1]-1,x,x[length(x)]+1)
polygon(x,chipInd,density=NA,col="#233142",lwd=2)
lines(x,predInd,col="#ff5959",lwd=2.5)
message("cat")
dev.off()
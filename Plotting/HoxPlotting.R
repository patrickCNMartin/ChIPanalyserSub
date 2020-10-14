#!/usr/bin/Rscript3.5.0

### Loading Libraries and Scripts

direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
library(ChIPanalyser)


#

###

# getting files
files<-dir()
flags<-c("HoX_ubx","HoX_abdb","HoX_dfd")
    # splitting files into their cats
    chips<-files[grepl("ChIP", files) & grepl(".Rda", files)]
    optimals<-files[grepl("Output",files) & grepl(".Rda", files)]

    chipselect<-c()
    optimalsselect<-c()
    for(i in seq_along(flags)){
        chipselect<-c(chipselect,chips[grep(flags[i],chips)])
        optimalsselect<-c(optimalsselect,optimals[grep(flags[i],optimals)])
    }


## selecting files with correct step size.
## you messed up the annotation of files
Chipstep100<-grep("step100", chipselect, value=T)
Optimalstep100<-grep("step100", optimalsselect, value=T)
### over validation
Chipstep100<-grep("Validation", chipselect, value=T)
Chipstep100<-grep("MSE", chipselect, value=T)
Optimalstep100<-grep("Validation", optimalsselect, value=T)
Optimalstep100<-grep("AUC", optimalsselect, value=T)


## next we split by tf

bytfChIP<-vector("list", length(flags))
bytfOptimal<-vector("list", length(flags))
names(bytfChIP)<-flags
names(bytfOptimal)<-flags

for(i in seq_along(bytfChIP)){
    bytfChIP[[i]]<-grep(flags[i],Chipstep100, value=T)
    bytfOptimal[[i]]<-grep(flags[i],Optimalstep100, value=T)

}

### extracting max values

training <- vector("list",length(bytfOptimal))
names(training)<-flags

for(i in seq_along(bytfOptimal)){
    training[[i]]<-vector("list",4)
    training[[i]][[1]]<-rep(0,length(bytfOptimal[[i]]))
    training[[i]][[2]]<-rep(0,length(bytfOptimal[[i]]))
    training[[i]][[3]]<-rep(0,length(bytfOptimal[[i]]))
    training[[i]][[4]]<-rep(0,length(bytfOptimal[[i]]))
    for(j in seq_along(bytfOptimal[[i]])){
        buffer<-get(load(bytfOptimal[[i]][[j]]))
        training[[i]][[1]][j]<-max(buffer[[1]][[2]][["AUC"]])
        training[[i]][[2]][j]<-max(buffer[[1]][[2]][["spearman"]])
        training[[i]][[3]][j]<-min(buffer[[1]][[2]][["MSE"]])
        training[[i]][[4]][j]<-max(buffer[[1]][[2]][["recall"]])
    }
    names(training[[i]])<-c("AUC","spearman","MSE","recall")
}


## extracting top metric over validation
validation <- vector("list",length(bytfOptimal))
names(validation)<-flags

for(i in seq_along(bytfOptimal)){
    validation[[i]]<-vector("list",1)
    validation[[i]][[1]]<-rep(0,length(bytfOptimal[[i]]))
    validation[[i]][[2]]<-rep(0,length(bytfOptimal[[i]]))
    validation[[i]][[3]]<-rep(0,length(bytfOptimal[[i]]))
    validation[[i]][[4]]<-rep(0,length(bytfOptimal[[i]]))
    for(j in seq_along(bytfOptimal[[i]])){
        buffer<-get(load(bytfOptimal[[i]][[j]]))
        loc1<-profiles(buffer$Gof)[[1]]
        validation[[i]][[1]][j] <- mean(sapply(loc1,"[","AUC")[10:60])
        validation[[i]][[2]][j]<-mean(sapply(loc1,"[","spearman")[10:60])
        validation[[i]][[3]][j]<-mean(sapply(loc1,"[","MSE")[10:60])
        validation[[i]][[4]][j]<-mean(sapply(loc1,"[","recall")[10:60])
    }
    names(validation[[i]])<-c("AUC","spearman","MSE","recall")
}


## posititon of best performing one

tops<-vector("list", length(training))

for(i in seq_along(tops)){
     maxiLoc<-which(training[[i]][[1]]==max(training[[i]][[1]]))
    tops[[i]]<-c(training[[i]][[1]][maxiLoc], bytfOptimal[[i]][maxiLoc])
}

## best over validation
tops<-vector("list", length(validation))

for(i in seq_along(tops)){
     maxiLoc<-which(validation[[i]][[1]]==max(validation[[i]][[1]]))
    tops[[i]]<-c(validation[[i]][[1]][maxiLoc], bytfOptimal[[i]][maxiLoc])
}



## first stage plotting
cols<-c("#233142","#ff5959","#4f9da6","#ff5959")
line <-1:3
accessThresholds<-seq(0,1,by=0.25)

plot(0,type="n", xlab='',ylab='',xlim=c(1,25),ylim=c(0.5,1),axes=F)
par(family="sans")
axis(1,at=seq(1,25,length.out=length(accessThresholds)),labels=accessThresholds,las=2)
axis(2, at=seq(0.5,1,by =0.1), labels=seq(0.5,1,by =0.1), las=2)
legend("topleft",legend=c("Ubx","Abd-b","Dfd"),col=cols,lty=1:3,bty="n",lwd=2,cex=1.4)
title(xlab="QDA")
title(ylab="AUC")


 for(i in seq_along(validation)){
    lines(seq_along(validation[[i]]$AUC),validation[[i]]$AUC,col=cols[i],lty=line[i] ,lwd=2)
 }


## making bars
extractAUC<-as.numeric(sapply(tops,"[[",1))
extractThreshold <- as.numeric(sapply(strsplit(sapply(tops,"[[",2), "_"),"[[",4))
names(extractThreshold)<-c("Ubx","Abd-b","Dfd")

par(family="sans")
par(xpd=NA)
par(mar=c(4,6,4,4))
barplot(extractThreshold, ylim=c(0,1),col=cols,cex.names=1.4,axes=F)
title(ylab="QDA Threshold",line=4.5, cex.lab=1.4)
axis(2,at=seq(0,1,by=0.1), seq(0,1,by=0.1),las=2,cex.axis=1.4)
title(ylab="QDA", cex.lab=2,line=6)
text(x=0.7,y=1.02, label=paste("AUC",round(extractAUC[1],digits=3)),cex=1.2)
text(x=1.9,y=0.78, label=paste("AUC",round(extractAUC[2],digits=3)),cex=1.2)
text(x=3.1,y=0.88, label=paste("AUC",round(extractAUC[3],digits=3)),cex=1.2)


################################################################################
################################################################################

flags<-c("HoX_ubx","HoX_abdb","HoX_dfd")
    # splitting files into their cats
    chips<-files[grepl("ChIP", files) & grepl(".Rda", files)]
    optimals<-files[grepl("Output",files) & grepl(".Rda", files)]

    chipselect<-c()
    optimalsselect<-c()
    for(i in seq_along(flags)){
        chipselect<-c(chipselect,chips[grep(flags[i],chips)])
        optimalsselect<-c(optimalsselect,optimals[grep(flags[i],optimals)])
    }



ChipValMSE<-chipselect[grepl("Validation", chipselect) & grepl("AUC",chipselect)]
OptimalValMSE<-optimalsselect[grepl("Validation", optimalsselect) & grepl("AUC",optimalsselect)]


## next we split by tf

bytfChIP<-vector("list", length(flags))
bytfOptimal<-vector("list", length(flags))
names(bytfChIP)<-flags
names(bytfOptimal)<-flags

for(i in seq_along(bytfChIP)){
    bytfChIP[[i]]<-grep(flags[i],ChipValMSE, value=T)
    bytfOptimal[[i]]<-grep(flags[i],OptimalValMSE, value=T)

}

selectionListOptimal<-c("HoX_ubx_ATACAccess_0.99_AUCOptimalOutputValidation.Rda",
                        "HoX_dfd_ATACAccess_0.7_AUCOptimalOutputValidation.Rda",
                      "HoX_abdb_ATACAccess_0.8_AUCOptimalOutputValidation.Rda")

selectionListChIP<-c("HoX_ubx_ATACAccess_0.99_AUCChIPValidation.Rda",
                        "HoX_dfd_ATACAccess_0.7_AUCChIPValidation.Rda",
                        "HoX_abdb_ATACAccess_0.8_AUCChIPValidation.Rda")
selectionAccess<-c("../objects/ATACAccess_0.99.Rda","../objects/ATACAccess_0.7.Rda","../objects/ATACAccess_0.8.Rda")

#for(i in seq_along(selectionListOptimal)){
  #  pdf(paste0(flags[i],".pdf"),width=14,height=4)
    #chipLoc<-get(load(selectionListChIP[i]))
  #  optimalLoc<-get(load(selectionListOptimal[i]))
  #  access<-get(load(selectionAccess[i]))
  #  plotOccupancyProfile(optimalLoc$ChIPProfile,chipLoc,chromatinState=access)
  #  dev.off()
#}






pdf("HoX_Profiles.pdf", width=10,height=12)
layout(matrix(c(1,2,3,3,4,4,5,5), ncol=2, byrow=T))
## first stage plotting
cols<-c("#233142","#ff5959","#4f9da6","#ff5959")
line <-1:3
accessThresholds<-seq(0,1,by=0.25)

plot(0,type="n", xlab='',ylab='',xlim=c(1,25),ylim=c(0.5,1),axes=F)
par(family="sans",xpd=NA)
axis(1,at=seq(1,25,length.out=length(accessThresholds)),labels=accessThresholds,las=2)
axis(2, at=seq(0.5,1,by =0.1), labels=seq(0.5,1,by =0.1), las=2)
legend(x=0,y=1.1,legend=c("Ubx","Abd-b","Dfd"),col=cols,lty=1:3,bty="n",lwd=2,cex=1.4)
text(x=-2,y=1.1, label="A",cex=4)
title(xlab="QDA")
title(ylab="AUC")

for(i in seq_along(training)){
   lines(seq_along(training[[i]]$AUC),training[[i]]$AUC,col=cols[i],lty=line[i] ,lwd=2)
}

## making bars
#extractAUC<-as.numeric(sapply(tops,"[[",1))
#extractThreshold <- as.numeric(sapply(strsplit(sapply(tops,"[[",2), "_"),"[[",4))
#names(extractThreshold)<-c("Ubx","Abd-b","Dfd")

par(family="sans")
par(xpd=NA)
par(mar=c(4,6,4,4))
barplot(extractThreshold, ylim=c(0,1),col=cols,cex.names=1.4,axes=F)
title(ylab="QDA Threshold",line=4.5, cex.lab=1.2)
axis(2,at=seq(0,1,by=0.1), seq(0,1,by=0.1),las=2,cex.axis=1.2)
#title(ylab="QDA", cex.lab=2,line=6)
text(x=0.7,y=1.04, label=paste("AUC",round(extractAUC[1],digits=3)),cex=1.2)
text(x=1.9,y=0.68, label=paste("AUC",round(extractAUC[2],digits=3)),cex=1.2)
text(x=3.1,y=0.78, label=paste("AUC",round(extractAUC[3],digits=3)),cex=1.2)
text(x=-0.35,y=1.15, label="B",cex=4)
labs <-LETTERS[3:6]
reg <- c(45,26,17)
tits <-c("Ubx - QDA 0.95 - Validation Set","Dfd - QDA 0.70 - Validation Set","Abd-b - QDA 0.80 - Validation Set")
for(j in seq_along(selectionListOptimal)){
  chipLoc<-get(load(selectionListChIP[j]))
  optimalLoc<-get(load(selectionListOptimal[j]))
  access<-get(load(selectionAccess[j]))

  setSequence<-loci(chipLoc)[reg[j]]
  ChIP<-scores(chipLoc)[[reg[j]]]
  pred<-profiles(optimalLoc$ChIPProfile)[[1]][[reg[j]]]

  cols<-c("#facf5a","#233142","#ff5959","#4f9da6")
  par(mar=c(5,2,4,1))
  par(family="sans")
  x<-seq(start(setSequence),end(setSequence),by=10)
  x<-c(x[1]-1,x,x[length(x)]+1)
  plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
  axis(1,at=seq(start(setSequence),end(setSequence),by=1000), labels=seq(start(setSequence),end(setSequence),by=1000))
  #axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)

  title(xlab=paste0("Genomic Position on : ",seqnames(setSequence)))
  title(main=tits[j])
  text(x=start(setSequence),y=1,labs[j], cex=4)

  ## plotting stuff
  noaccess<-.AccessExtract(setSequence,access)[[1]]
  for(i in seq_len(nrow(noaccess))){
      rect(noaccess[i,"start"],0,noaccess[i,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
      #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
  }

  chipInd<-c(0,ChIP[seq(0,length(ChIP),by=10)],0)
  predInd<-c(0,pred$ChIP,0)

  polygon(x,chipInd,col=cols[2],lwd=2)
  lines(x,predInd,col=cols[3],lwd=2)



}

dev.off()












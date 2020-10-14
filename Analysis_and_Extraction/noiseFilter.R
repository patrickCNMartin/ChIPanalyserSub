###############################################################################
                        #Noise filter plotting #
###############################################################################

## Data set up starting with the usual

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


## original colours for plotting
#cols<-c("#4f9da6","#233142","#ff5959","#facf5a")
cols2<-c("#233142","#4f9da6","#ff5959","#facf5a")

colfunc<-colorRampPalette(c("#4f9da6","#233142"))
cols<-colfunc(4)
## first let's load the data
## Lets start with noisy data
zero<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005zero_ChIPTraining.Rda"))
mean<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005mean_ChIPTraining.Rda"))
median<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005median_ChIPTraining.Rda"))
sigmoid<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005sigmoid_ChIPTraining.Rda"))

noisy<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
noisy<-lapply(noisy, scores)
noisy<-lapply(noisy,function(x,y){x[y]},c(8,10))

#lets get clean data
zero<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005zero_ChIPTraining.Rda"))
mean<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005mean_ChIPTraining.Rda"))
median<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005median_ChIPTraining.Rda"))
sigmoid<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005sigmoid_ChIPTraining.Rda"))

clean<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
clean<-lapply(clean, scores)
clean<-lapply(clean,function(x,y){x[y]},c(4,5))
# lets get combi Data
zero<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce10_new_S2_DHS_005zero_ChIPTraining.Rda"))
mean<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce10_new_S2_DHS_005mean_ChIPTraining.Rda"))
median<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce10_new_S2_DHS_005median_ChIPTraining.Rda"))
sigmoid<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce10_new_S2_DHS_005sigmoid_ChIPTraining.Rda"))

combi<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
combi<-lapply(combi,scores)
combi<-lapply(combi,function(x,y){x[y]},c(7,2))
## lets sub the shit out of these methods_ChIPanalyser

chip<-list("Noisy"=noisy,"Clean"=clean,"Combined"=combi)
## Creating setSequence

generateSetSequence<-function(chip){
    name<-names(chip[[1]])

    chr<-sapply(strsplit(name,":"),"[[",1)
    start<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",1)
    end<-sapply(strsplit(sapply(strsplit(name,":"),"[[",2),"\\.."),"[[",2)
    setSequence<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))
    return(setSequence)
}

setSequence<-lapply(chip,generateSetSequence)

## lets do some stats
optimalParam<-function(noiseFilters,method="AUCMean",method1="AUC"){

   # extracting top paramters to get the adequate parameters based on method
   topParam<-as.numeric(noiseFilters[[1]][[1]][[method]])

   # let's get accuracy score for top paramters
   accu <- searchSites(noiseFilters[[4]],topParam[1],topParam[2])
   accu <- sapply(accu[[1]],function(x){return(x[[1]][method1])})
   return(list(topParam,accu))
}
# load that optimal Data
## Lets start with noisy data
zero<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005zero_OptimalOutputTraining.Rda"))
mean<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005mean_OptimalOutputTraining.Rda"))
median<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005median_OptimalOutputTraining.Rda"))
sigmoid<-get(load("singleNoisy/BG3_modEncode_3674_CTCF_AllLoci_reduce10_new_BG3_DHS_005sigmoid_OptimalOutputTraining.Rda"))

noisy<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
noisy<-lapply(noisy,function(x){buf<-profiles(x$goodnessOfFit)[[1]];sapply(buf, "[[","AUC")})


#lets get clean data
zero<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005zero_OptimalOutputTraining.Rda"))
mean<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005mean_OptimalOutputTraining.Rda"))
median<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005median_OptimalOutputTraining.Rda"))
sigmoid<-get(load("singleClean/S2_modEncode_2639_CTCF_AllLoci_reduce10_new_S2_DHS_005sigmoid_OptimalOutputTraining.Rda"))

clean<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
clean<-lapply(clean,function(x){buf<-profiles(x$goodnessOfFit)[[1]];sapply(buf, "[[","AUC")})


# lets get combi Data
zero<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce3293_S2_DHS_005zero_OptimalOutputTraining.Rda"))
mean<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce3293_S2_DHS_005mean_OptimalOutputTraining.Rda"))
median<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce3293_S2_DHS_005median_OptimalOutputTraining.Rda"))
sigmoid<-get(load("ctcfCombi/S2_CTCFcombi_AllLoci_reduce3293_S2_DHS_005sigmoid_OptimalOutputTraining.Rda"))

combi<-list("Zero"=zero,"Mean"=mean,"Median"=median,"Sigmoid"=sigmoid)
combi<-lapply(combi,function(x){buf<-profiles(x$goodnessOfFit)[[1]];sapply(buf, "[[","AUC")})


Noisefilter<-list(noisy,clean,combi)

## First set of plots

pdf("NoiseFilter_all.pdf", height=12,width=20)
layout(matrix(1:9,ncol=3,byrow=T),width=c(9,9,4))
par(oma=c(0,0,3,0))

par(xpd=NA)
par(family="sans")
biglabs<-LETTERS[1:6]
count<-1
for(i in seq_along(chip)){

   for(j in seq_along(setSequence[[1]])){
      x<-seq(start(setSequence[[i]])[j],end(setSequence[[i]])[j],by=10)
      x<-c(x[1]-1,x,x[length(x)]+1)
      par(xpd=NA)
      if(j==1){
         par(mar=c(5,12,1,0))
      }else{
        par(mar=c(5,1,1,4))
      }
      plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(-3.3,1))
      axis(1,at=seq(start(setSequence[[i]])[j]+2000,end(setSequence[[i]])[j],by=2000), labels=seq(start(setSequence[[i]])[j]+2000,end(setSequence[[i]])[j],by=2000),cex.axis=1.2)

      if(j==1){
        text(x=start(setSequence[[i]])[j]-3000,y=1.1, labels=biglabs[count],cex=4)
        axis(2,at=c(-3.3,-2.2,-1.1,0),labels=c("Sigmoid","Median","Mean","Zero"),las=2,cex.axis=2)
        title(ylab=names(chip)[i],line=9, cex.lab=2,font=2)
      }
      title(xlab=paste0("Genomic Position on : ",seqnames(setSequence[[i]][j])),font=4,cex.lab=1.6)


      zeorInd<-c(0,chip[[i]][[1]][[j]][seq(0,length(chip[[i]][[1]][[j]]),by=10)],0)
      meanInd<-c(0,chip[[i]][[2]][[j]][seq(0,length(chip[[i]][[2]][[j]]),by=10)],0)
      medianInd<-c(0,chip[[i]][[3]][[j]][seq(0,length(chip[[i]][[3]][[j]]),by=10)],0)
      sigmoidInd<-c(0,chip[[i]][[4]][[j]][seq(0,length(chip[[i]][[4]][[j]]),by=10)],0)


      polygon(x,zeorInd,density=NA,col=cols[1],lwd=2)
      polygon(x,meanInd-1.1, density=NA,col=cols[2],lwd=2)
      polygon(x,medianInd-2.2, density=NA,col=cols[3],lwd=2)
      polygon(x,sigmoidInd-3.3, density=NA,col=cols[4],lwd=2)
    }
    par(xpd=NA)
    par(mar=c(7,4,6,2))
    boxplot(Noisefilter[[i]],col=cols2,frame=F,main="AUC Distribution",cex.main=1.6,cex.axis=1.6,las=2,ylim=c(0,1))
    count<-count+1
    #print(max(Noisefilter[[i]][[4]]))
    if(i==1){
        #text(y=(max(Noisefilter[[i]][[4]],na.rm=T)+(max(Noisefilter[[i]][[4]],na.rm=T)-min(Noisefilter[[i]][[4]],na.rm=T))/4),x=(-0.3),biglabs[count],cex=4,xpd=NA)
        text(y=1.3,x=(-0.3),biglabs[count],cex=4,xpd=NA)
    } else if(i==2){
        #text(y=(max(Noisefilter[[i]][[4]],na.rm=T)+(max(Noisefilter[[i]][[4]],na.rm=T)-min(Noisefilter[[i]][[4]],na.rm=T))/1.5),x=(-0.3),biglabs[count],cex=4,xpd=NA)
        text(y=1.3,x=(-0.3),biglabs[count],cex=4,xpd=NA)
    } else {
    #text(y=(max(Noisefilter[[i]][[4]],na.rm=T)+(max(Noisefilter[[i]][[4]],na.rm=T)-min(Noisefilter[[i]][[4]],na.rm=T))),x=(-0.3),biglabs[count],cex=4,xpd=NA)
    text(y=1.3,x=(-0.3),biglabs[count],cex=4,xpd=NA)
    }
    count<-count+1
}
#mtext("Noise Filtering Methods - CTCF - Combi",outer=T,cex=2)

dev.off()

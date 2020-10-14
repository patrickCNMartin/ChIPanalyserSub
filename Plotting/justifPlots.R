###############################################################################
                          # Justification plots #
###############################################################################

## start off with usual loading

# Loading libraris and sourcing code
direc <- getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)


## sourcing scripts for analysis
setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyser_1.1")
files <- dir()
for (i in files) source(i)
setwd(direc)


## Loading data

## Loading DNA access
Access<-get(load("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/Data/DHS_500bp_BG3.Rda"))



## Lets got for some BG3 and BEAF
beaf32_3665_ChIP <- get(load("peakWindowreduce20/BG3_modEncode_3665_BEAF-32_reduce100sigmoid_AUCMeanChIP.Rda"))
beaf32_3665_Optimal <- get(load("peakWindowreduce20/BG3_modEncode_3665_BEAF-32_reduce100sigmoid_AUCMeanoptimalOutput.Rda"))

## lets also go for someBG3 ctcfCombi
ctcf_3280_ChIP <- get(load("peakWindowreduce20/BG3_modEncode_3280_CTCF_reduce100sigmoid_AUCMeanChIP.Rda"))
ctcf_3280_Optimal <- get(load("peakWindowreduce20/BG3_modEncode_3280_CTCF_reduce100sigmoid_AUCMeanoptimalOutput.Rda"))

## abdb original data
abdb_Optimal<-get(load("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/peakWindowreduce100/Kc167_abdb_reduce100optimalOutput.Rda"))
opti_mat<-abdb_Optimal[[2]]


## Let's get the profiles of interest

beafLoci<-c("chrX:17780001..17800000","chr3R:24860001..24880000")

ctcfLoci<-c("chr2R:18020001..18040000","chr2L:5260001..5280000")

beaf32_3665_ChIP <- beaf32_3665_ChIP[beafLoci]
ctcf_3280_ChIP <- ctcf_3280_ChIP[ctcfLoci]

optimalBeaf <- as.numeric(beaf32_3665_Optimal[[1]][[1]][["pearsonMean"]])
optimalCtcf <- as.numeric(ctcf_3280_Optimal[[1]][[1]][["pearsonMean"]])

chipBeaf<-searchSites(beaf32_3665_Optimal[[3]],optimalBeaf[1],optimalBeaf[2])[[1]][beafLoci]
chipBeaf<-lapply(chipBeaf, function(x){x$ChIP})
chipctcf<-searchSites(ctcf_3280_Optimal[[3]],optimalCtcf[1],optimalCtcf[2])[[1]][ctcfLoci]
chipctcf<-lapply(chipctcf, function(x){x$ChIP})

accuBeaf<-searchSites(beaf32_3665_Optimal[[4]],optimalBeaf[1],optimalBeaf[2])[[1]][beafLoci]
accuBeaf<-sapply(accuBeaf, function(x){x[[1]]["pearson"]})
accuctcf<-searchSites(ctcf_3280_Optimal[[4]],optimalCtcf[1],optimalCtcf[2])[[1]][ctcfLoci]
accuctcf<-sapply(accuctcf, function(x){x[[1]]["pearson"]})

## setting up setSequence
nameBEAF<-names(beaf32_3665_ChIP)
chr<-sapply(strsplit(nameBEAF,":"),"[[",1)
start<-sapply(strsplit(sapply(strsplit(nameBEAF,":"),"[[",2),"\\.."),"[[",1)
end<-sapply(strsplit(sapply(strsplit(nameBEAF,":"),"[[",2),"\\.."),"[[",2)
setSequenceBEAF<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))

nameCTCF<-names(ctcf_3280_ChIP)
chr<-sapply(strsplit(nameCTCF,":"),"[[",1)
start<-sapply(strsplit(sapply(strsplit(nameCTCF,":"),"[[",2),"\\.."),"[[",1)
end<-sapply(strsplit(sapply(strsplit(nameCTCF,":"),"[[",2),"\\.."),"[[",2)
setSequenceCTCF<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))
names(setSequenceCTCF)<-nameCTCF


setSequence<-c(setSequenceBEAF,setSequenceCTCF)

ChIP<-c(beaf32_3665_ChIP,ctcf_3280_ChIP)
chip<-c(chipBeaf,chipctcf)
accu<-c(accuBeaf,accuctcf)
tits<-c(rep("BEAF-32 - BG3",2),rep("CTCF - BG3",2))

NoAccess <- .AccessExtract(setSequence, Access)
### fuck it lets do the plots


cols<-c("#facf5a","#233142","#ff5959")
#cols2<-c("#233142","#4f9da6","#ff5959","#facf5a")

pdf("Justification.pdf",height=14,width=22)
layout(matrix(c(1,1,1,1,1,1,3,3,3,3,3,3,2,2,2,2,2,2,4,4,4,4,4,4,5,5,5,6,7,7,7,8,9,9,9,10), ncol=12, byrow=T),height=c(1.8,1.8,3.5))
par(xpd=NA)
biglabs<-c("A","B")
count<-1
for(i in seq_along(ChIP)){

  par(mar=c(6,4,4,1))
  par(family="mono")
  x<-seq(start(setSequence)[i],end(setSequence)[i],by=10)
  x<-c(x[1]-1,x,x[length(x)]+1)
  plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
  axis(1,at=seq(start(setSequence)[i]+2000,end(setSequence)[i],by=2000), labels=seq(start(setSequence)[i]+2000,end(setSequence)[i],by=2000),cex.axis=1.6)
  title(main=tits[i], cex.main=2,line=2)

  title(xlab=paste0("Genomic Position on : ",seqnames(setSequence[i])),cex.lab=1.6)
  if(i%%2!=0){
     text(x=start(setSequence)[i]-800,y=1,labels=biglabs[count],cex=5)
     count<-count+1
  }

  ## plotting stuff
  noaccess<-NoAccess[[i]]
  for(j in seq_len(nrow(noaccess))){
      rect(noaccess[j,"start"],0,noaccess[j,"end"],0.8,col=cols[1],density=10,angle=45,lwd=0.8,border=NA)
      rect(noaccess[j,"start"],0,noaccess[j,"end"],0.8,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
      rect(noaccess[j,"start"],0,noaccess[j,"end"],0.8,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
  }

  chipInd<-c(0,ChIP[[i]][seq(0,length(ChIP[[i]]),by=10)],0)
  predInd<-c(0,chip[[i]],0)

  polygon(x,chipInd,density=NA,col=cols[2],lwd=2)
  lines(x,predInd,col=cols[3],lwd=2.5)
  text(x=end(setSequence)[i]-3900,y=0.85,labels=paste("Pearson Correlation",round(accu[i],digits=3)),cex=1.6,font=2)


}

xlabs<-as.numeric(colnames(opti_mat[[1]]))
ylabs<-as.numeric(rownames(opti_mat[[1]]))
for(i in seq_along(opti_mat)){


        colfunc<-colorRampPalette(c("white",cols[i]))


        Colors <- colfunc(20)
        legend_image<-as.raster(matrix(rev(Colors),ncol=1))
        par(mar=c(7.5,7,4.5,1)+0.1)

        graphics::image(1:length(xlabs),1:length(ylabs),t(opti_mat[[i]]),
            axes = FALSE, xlab=" ", ylab=" ",col=Colors)
        box()

        if(i ==1){
            text(x=(-1),y=nrow(opti_mat[[i]])+2,labels="C", cex=5)
        }

        title(main=names(opti_mat)[i],cex.main=2)
        title( ylab="Scaling Factor", line=5, cex.lab=1.8)
        title(xlab="Number of Bound Molecules",line=5.5, cex.lab=1.8)
        axis(1,at=seq_along(xlabs),labels=F,cex.axis=1.4)
        text(seq_along(xlabs),y =(-0.2), srt = 45, adj = 1,labels = xlabs, xpd = TRUE,cex=1.4)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,las= HORIZONTAL<-1,cex.axis=1.4)
        #text(x=-2,y=length(ylabs)+2,labels=right[i],cex=4)
        # raster scacle
        par(mar=c(5.5,0,3.5,0.5))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
        text(x=1.6, y =seq(0,1,l=5) , labels = round(seq(min(opti_mat[[i]]),max(opti_mat[[i]]),l=5),2),cex=1.6)
        rasterImage(legend_image, 0, 0, 1,1)



}




dev.off()

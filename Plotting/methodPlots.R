

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

colsheat<-rep(c("#4f9da6","#233142","#ff5959"),times=2)
cols<-c(rep("#4f9da6",times=4),rep("#233142",times=8))

mains<-c("CTCF - Clean (modEncode 2639)","CTCF - Noisy (modEncode 3674)","CTCF - Combined")
pdf("methods_ChIPanalyserHeat.pdf",width=12,height=16)

#layout(matrix(cbind(c(1,2,3),c(4,6,8),c(5,7,9)),ncol=3),width=c(8,3,0.8),height=c(1,1,1))
#layout(matrix(cbind(c(1,3,5),c(2,4,6),c(7,8,9)),ncol=3),width=c(3,0.8,8),height=c(1,1,1))
layout(matrix(cbind(c(1,3,5),c(2,4,6),c(7,9,11),c(8,10,12)),ncol=4),width=c(3,0.8,3,0.8),height=c(1,1,1))

par(xpd=NA)
par(oma=c(0,0,0,3))
par(family="sans")
ylabs<-as.numeric(rownames(heats[[1]]))
xlabs<-as.numeric(colnames(heats[[1]]))
right<-LETTERS[1:6]
for(i in seq_along(heats)){
        #if(grepl("MSE",mainTitle[i])|grepl("geometric",mainTitle[i])|grepl("ks", mainTitle[i])){
            #colfunc<-colorRampPalette(c(cols[i],"white"))
        #}else{
            colfunc<-colorRampPalette(c("white",colsheat[i]))
      #  }

        Colors <- colfunc(20)
        legend_image<-as.raster(matrix(rev(Colors),ncol=1))
        if(i<=3){par(mar=c(5.5,7,4.5,1)+0.1)}else{
        par(mar=c(5.5,6,4.5,1)+0.1)
        }
        graphics::image(1:length(xlabs),1:length(ylabs),t(heats[[i]]),
            axes = FALSE, xlab=" ", ylab=" ",col=Colors)
        box()
        title(main=names(heats)[i],cex.main=2)
        if(i <= 3){title( ylab="Scaling Factor", line=4.5, cex.lab=1.4)}
        title(xlab="Number of Bound Molecules",line=4.5, cex.lab=1.4)
        axis(1,at=seq_along(xlabs),labels=F,cex.axis=1.4)
        text(seq_along(xlabs),y =(-0.2), srt = 45, adj = 1,labels = xlabs, xpd = TRUE,cex=1.4)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,las= HORIZONTAL<-1,cex.axis=1.4)
        text(x=-2,y=length(ylabs)+2,labels=right[i],cex=4)
        # raster scacle
        par(mar=c(5.5,0,3.5,0.5))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
        text(x=1.6, y =seq(0,1,l=5) , labels = round(seq(min(heats[[i]]),max(heats[[i]]),l=5),2),cex=1.4)
        rasterImage(legend_image, 0, 0, 1,1)
    }


dev.off()
colsheat<-rep(c("#4f9da6","#233142","#ff5959"),times=2)
cols<-c(rep("#4f9da6",times=4),rep("#233142",times=8))

mains<-c("CTCF - Clean (modEncode 2639)","CTCF - Noisy (modEncode 3674)","CTCF - Combined")
pdf("methods_ChIPanalyserProfile.pdf",width=12,height=16)
par(mfrow=c(3,1))
par(xpd=NA)
par(oma=c(0,0,0,3))
par(family="sans")
left<-LETTERS[1:3]
for( j in seq_along(ChIP)){
   x<-seq(start(setSequence)[j],end(setSequence)[j],by=100)
   x<-c(x[1]-1,x,x[length(x)]+1)

   name<-c("ChIP-seq","MSE","K-S Distance","Geometric","Recall","Pearson","Spearman","Kendall","Precision","F-score","Accuracy","MCC","AUC")
   lims<-seq(0,by=(-1.2),length.out=length(name))
   par(mar=c(5,10.5,4,9.5))
   par(xpd=T)
   plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(min(lims),1))

   axis(1,at=seq(start(setSequence)[j],end(setSequence)[j],by=1000), labels=seq(start(setSequence)[j],end(setSequence)[j],by=1000),cex.axis=1.4)

   axis(2,at=lims,labels=name,las=2,cex.axis=1.7)


   title(xlab=paste0("Genomic Position on : ",seqnames(setSequence[j])),cex.lab=1.4,font=2)
   title(main=mains[j],cex.main=2)
   text(x=x[1]-3000,y=2.6,labels=left[j],cex=4)
   #mtext(side=4
   #text(
   chipsInd<-c(0,ChIP[[j]][seq(0,length(ChIP[[j]]),by=100)],0)
   polygon(x,chipsInd,density=NA,col="#ff5959",lwd=2)
   pos<-c(5.5,5,6,8,1,2,3,7,9,10,11,12)
   buffername<-c("MSE","KsDist","geometric","recall","pearson","spearman","kendall","precision","f1","accuracy","MCC","AUC")

   estiInternal<-estimates[[j]]
   predictionsInternal<-predictions[[j]]
   optimalInternal<-optimal[[j]]

   qual<-rep(0,length(buffername))
   limint<-lims[2:length(lims)]
   for(i in seq_along(buffername)){
       #print(names(optimalInternal[[1]]))
       paraset<-optimalInternal[[1]][[buffername[i]]]
       if(class(paraset)=="matrix"){

           paraset<-c(paraset[nrow(paraset),1],paraset[nrow(paraset),1])
       }
       #if(buffername[i]=="AUC")paraset<-c(0.5,5000)

       paraset<-as.numeric(paraset)
       #predicted<-searchSites(predictionsInternal,lambdaPWM=paraset[1],BoundMolecules=paraset[2],names(ChIP)[j])
       predicted<-searchSites(predictionsInternal,lambdaPWM=paraset[1],BoundMolecules=paraset[2],names(ChIP)[j])
       predicted<-c(0,predicted[[1]][[1]]$ChIP,0)

       polygon(x,predicted+limint[i], density=NA,col=cols[i],lwd=2)

       buffer<-searchSites(estiInternal,lambdaPWM=paraset[1],BoundMolecules=paraset[2],names(ChIP)[j])
       if(!is.na(buffer[[1]][[1]][buffername[i]])){
          qual[i]<-buffer[[1]][[1]][buffername[i]]

       } else{
          qual[i]<-NaN
       }

   }

   axis(4,at=limint,labels=signif(qual,digits=4),las=2,cex.axis=1.7,line=-1)

   #mtext("Quality Assessment Method - CTCF - Clean",outer=T,cex=1.8,font=2)
   mtext("Associated Score",side=4,outer=T ,cex=1.4,font=2)

}

dev.off()

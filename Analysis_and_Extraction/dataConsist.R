
# Data set up starting with the usual

# Loading libraris and sourcing code

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)

library(ChIPanalyser)

source("DataHand.R")



### based on script in dataLoad
# if you want to redo this you have to have the outcome of
##  data <- bytf(nameMatch(AUCdhs,AUCnull))
## retuns matrix with ordered TF and duplicate drop, unwanted drop


overlayMat<-function(mat, method="geometricMean"){
   over<-matrix(0,nrow=14,ncol=17)
   mats<-vector("list",ncol(mat))
   for(i in seq_along(mats)){
     buffer<-mat[,i]
     dim(buffer)<-c(14,17)
     mats[[i]]<-buffer
   }

   for(k in seq_along(mats)){
     buffer<-mats[[k]]
     idx<-seq_along(as.vector(mats[[k]]))


     # define top hits
     top<-head(idx,floor(length(idx)/5))

     # Top hits
     if(method %in% c("geometricMean","MSEMean","ksMean","recallMean")){
         ord<-order(buffer,decreasing=F)
     } else{
         ord<-order(buffer,decreasing=T)
     }
     ##creating a ordered value matrix
     optimalMatrix<-matrix(match(idx,ord),ncol=17, nrow=14)

     ## Location of top hits in ordered value matrix
     topHits<-which(optimalMatrix<max(top),arr.ind=T)

     # Creating empty matrix
     mat<- matrix(0,ncol=17,nrow=14)

     ## Keep in mind that there are more Similarity methods
     ## This will pull it towards those optimal Parameters
     ## For better results balence it
     ## or only take two Parameters

     ## Overlaying top hits
     mat[topHits]<-1


     ## Adding top hits to matrix
     over<-over+mat
     }
     over<-over/max(over)
    colnames(over)<-c(1, 10, 20, 50, 100,
         200, 500,1000,2000, 5000,10000,20000,50000, 100000,
         200000, 500000, 1000000)
    rownames(over)<-c(0.25, 0.5, 0.75, 1, 1.25,
        1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5)
     return(over)

}

bytfmat<-function(mat,tf=c("CTCF","BEAF","Hw"),cell=c("BG3","Kc167","S2")){
   res<-vector("list",length(tf)*length(cell))
   count<-1
   tag<-c()
   for(i in seq_along(tf)){
      for(j in seq_along(cell)){
       res[[count]]<-mat[,which(grepl(tf[i],colnames(mat)) & grepl(cell[j],colnames(mat)))]
       tag<-c(tag,paste(tf[i],cell[j]))
       count<-count+1
     }
   }
   names(res)<-tag
   return(res)
}


decayMse<-bytf(nameMatch(MSEdhs,MSEnull))
decayAUC<-bytf(nameMatch(AUCdhs,AUCnull))
decaygeo<-bytf(nameMatch(GEOdhs,GEOnull))




submatsAUC<-bytfmat(decayAUC[[1]])
submatsGEO<-bytfmat(decaygeo[[1]])
submatsMSE<-bytfmat(decayMse[[1]])

overlayMatricesAUC<-lapply(submatsAUC,overlayMat,method="AUCMean")
overlayMatricesGEO<-lapply(submatsGEO,overlayMat,method="geometricMean")
overlayMatricesMSE<-lapply(submatsMSE,overlayMat,method="MSE")

## overlay of methods (only two to avoid biased towards certain type of methods)


overlayMatrices<-vector("list", length(overlayMatricesAUC))
for(i in seq_along(overlayMatrices)){
    overlayMatrices[[i]]<-overlayMatricesAUC[[i]]+overlayMatricesMSE[[i]]

}
names(overlayMatrices)<-names(overlayMatricesAUC)

overlayMatrices<-overlayMatricesAUC
overlayMatrices<-overlayMatricesGEO
overlayMatrices<-overlayMatricesMSE

names(overlayMatrices)<-gsub("Hw","su(Hw)",names(overlayMatrices))

colfunc<-colorRampPalette(c("white","#ff5959","#233142"))
labs<-LETTERS[1:9]



pdf("overlayMatriciesGEO.pdf",height=18,width=18)
par(family="mono")
par(xpd=NA)
layout(matrix(1:18,ncol=6,byrow=T),width=rep(c(1,0.28),times=3))
xlabs<-c(1, 10, 20, 50, 100,
    200, 500,1000,2000, 5000,10000,20000,50000, 100000,
    200000, 500000, 1000000)
ylabs<-c(0.25, 0.5, 0.75, 1, 1.25,
    1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5)

for(i in seq_along(overlayMatrices)){
par(mar=c(8,10,8,0))
Colors <- colfunc(15)
legend_image<-as.raster(matrix(rev(Colors),ncol=1))
graphics::image(1:length(xlabs),1:length(ylabs),t(overlayMatrices[[i]]),
axes = FALSE, xlab=" ", ylab=" ",col=Colors)
text(x=(-1), y=length(ylabs)+2,labels=labs[i],cex=5)
title(main=names(overlayMatrices)[i],cex.main=2)
title( ylab="Scaling Factor", line=6, cex.lab=1.6)
title(xlab="Number of Bound Molecules",line=6.5, cex.lab=1.6)
axis(1,at=seq_along(xlabs),labels=F)
text(seq_along(xlabs),y =(-0.2), srt = 45, adj = 1,labels = xlabs, xpd = TRUE,cex=1.4)
axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,las= HORIZONTAL<-1,cex.axis=1.4)
par(mar=c(5,1,6,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.6, y =seq(0,1,l=5) , labels = round(seq(min(overlayMatrices[[i]]),max(overlayMatrices[[i]]),l=5),2),cex=1.4)
rasterImage(legend_image, 0, 0, 1,1)
}
dev.off()

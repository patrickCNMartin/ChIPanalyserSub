tf<-c("CTCF","BEAF-32","su(Hw)")
cell<-c("Kc167","BG3","BG3","S2","Kc167","S2")
labs<-LETTERS[1:6]
## you need to set up your data, which is not here because you are an idiot :D
pdf(paste0("RNA_rescale_ChIPanalyser_","MSE",".pdf"), width=27,height=15)
#par(oma=c(0,0,9,0))
layout(matrix(cbind(c(1,2,5,6,9,10),c(3,4,7,8,11,12),c(13,13,14,14,15,15)),ncol=3), width=c(7,7,5.5),height=c(1,1,1))
par(family="sans")
par(xpd=NA)

cols<-c("#ff5959","#233142","#facf5a","#facf5a")
count<-0

for(i in seq_along(rescaledExtract)){
    setTrain <- loci(ChIPProfiles[[i]])
    scoresTrain<- scores(ChIPProfiles[[i]])
print(count)
    param <- as.numeric(predictions[[i]][[1]][[1]][[method]])

    lambda <- param[1]
    bound <-param[2]

    count<-count+1
subcount<-1
  for(k in TFsub[[i]]){

    predictionTrain <- searchSites(predictions[[i]]$ChIPProfiles, lambda, bound,names(scoresTrain)[k])
    x<-seq(start(setTrain)[k],end(setTrain)[k],by=100)
    x<-c(x[1]-1,x,x[length(x)]+1)
    if(subcount==1){par(mar=c(5,2,4.5,2))}else{par(mar=c(5,2,3.8,2))}

    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(setTrain)[k],end(setTrain)[k]),ylim=c(0,1))

    title(xlab=paste0("Genomic Position on ",as.character(seqnames(setTrain))[k]),cex.lab=1.5)

    if(subcount==1){
    title(main=paste0(tf[i]," in ",cell[count]," - lambda = ",lambda," & Bound Molecules = ",bound),cex.main=1.8,line=0.5)
    text(x=(x[1]-1000), y=1.2, labels=labs[count],cex=4)
    }
    subcount<-subcount+1
    axis(1,at=round(seq(start(setTrain)[k],end(setTrain)[k],length.out=10)),labels=round(seq(start(setTrain)[k],end(setTrain)[k],length.out=10)),cex.axis=1.5)

    noaccess<-.AccessExtract(setTrain[k],AccessOriginal[[i]])[[1]]

    for(j in seq_len(nrow(noaccess))){
        rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
    }

     local<-scoresTrain[[k]]
    chipInd<-c(0,local[seq(0,length(local),by=100)],0)
    predInd<-c(0,predictionTrain[[1]][[1]]$ChIP,0)

    polygon(x,chipInd,density=NA,col="#233142",lwd=2)
    lines(x,predInd,col="#ff5959",lwd=2.5)
  }
  count<-count+1
  print(count)
  subcou<-1
  for(k in TFsub[[i]]){
  boundres <- as.numeric(round((param[2]/scale[[i]][[1]])*scale[[i]][[2]]))

    validationScore <- scores(rescaledChIPsignal[[i]])
    validationSet <-loci(rescaledChIPsignal[[i]])

    x<-seq(start(validationSet)[k],end(validationSet)[k],by=100)
    x<-c(x[1]-1,x,x[length(x)]+1)
    if(subcou==1){par(mar=c(5,2,4.5,2))}else{par(mar=c(5,2,3.8,2))}

    plot(0,type="n", axes=FALSE,xlab="",ylab="",xlim=c(start(validationSet)[k],end(validationSet)[k]),ylim=c(0,1))

    title(xlab=paste0("Genomic Position on ",as.character(seqnames(validationSet))[k]),cex.lab=1.5)

    if(subcou==1){
    title(main=paste0(tf[i]," in ",cell[count]," - lambda = ",lambda," & Bound Molecules = ",boundres),cex.main=1.8,line=0.5)
    text(x=(x[1]-1000), y=1.2, labels=labs[count],cex=4)
    }
    subcou<-subcou+1
    axis(1,at=round(seq(start(validationSet)[k],end(validationSet)[k],length.out=10)),labels=round(seq(start(validationSet)[k],end(validationSet)[k],length.out=10)),cex.axis=1.5)

    noaccess<-.AccessExtract(validationSet[k],Access[[i]])[[1]]

    for(j in seq_len(nrow(noaccess))){
        rect(noaccess[j,"start"],0,noaccess[j,"end"],0.9,col="#facf5a",density=50,angle=45,lwd=1,border=NA)
        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=135,lwd=1,border=NA)
        #rect(noaccess[[k]][j,"start"],0,noaccess[[k]][j,"end"],0.9,col="#facf5a",density=10,angle=90,lwd=1,border=NA)
    }

     local<-validationScore[[k]]
    chipInd<-c(0,local[seq(0,length(local),by=100)],0)
    Carry<-profiles(rescaledExtract[[i]][[1]]$ChIPPrediction)
    CarryInd<-c(0,Carry[[1]][[k]]$ChIP,0)

    rescale <-profiles(rescaledExtract[[i]][[2]]$ChIPPrediction)
    rescaleInd<-c(0,rescale[[1]][[k]]$ChIP,0)

    fold10 <-profiles(rescaledExtract[[i]][[3]]$ChIPPrediction)
    foldInd10<-c(0,fold10[[1]][[k]]$ChIP,0)

    fold100 <-profiles(rescaledExtract[[i]][[4]]$ChIPPrediction)
    foldInd100<-c(0,fold100[[1]][[k]]$ChIP,0)


    polygon(x,chipInd,density=NA,col="#233142",lwd=2)
    lines(x,CarryInd,col="#56B4E9",lwd=2.5)
    lines(x,rescaleInd,col="#ff5959",lwd=2.5,lty=2)
    lines(x,foldInd10,col="#a64f9d",lwd=2.5,lty=1)
    lines(x,foldInd100,col="#009E73",lwd=2.5,lty=1)
    legend(x=x[length(x)],y=1,col=c("#56B4E9","#ff5959","#a64f9d","#009E73"),
    lty=c(1,2,1),legend=c("Carry-Over","Rescaled","10 Fold","100 Fold"),bty="n",lwd=rep(2.5,3),cex=2)
  }


}






## initating boxplot for loop
labs<-LETTERS[7:9]
par(xpd=NA)
for(i in seq_along(rescaledExtract)){
    param<-as.numeric(predictions[[i]][[1]][[1]][[method]])
    ext <- searchSites(predictions[[i]]$goodnessOfFit,param[1],param[2])[[1]]
    ext <- sapply(ext, function(x){return(x[method])})

    carry <- profiles(rescaledExtract[[i]][["Carry"]][["GoF"]])[[1]]
    carry <- sapply(carry, function(x){return(x[method])})
    rescale <- profiles(rescaledExtract[[i]][["Rescaled"]][["GoF"]])[[1]]
    rescale <- sapply(rescale, function(x){return(x[method])})
    fold <- profiles(rescaledExtract[[i]][["100Fold"]][["GoF"]])[[1]]
    fold <- sapply(fold, function(x){return(x[method])})
    fold10 <- profiles(rescaledExtract[[i]][["10Fold"]][["GoF"]])[[1]]
    fold10 <- sapply(fold10, function(x){return(x[method])})
    par(xpd=NA)
    par(mar=c(4,17,4,2))
    dat<-list("Estimated"=ext,"Rescaled"=rescale,"Carry-Over"=carry,"10 Fold"=fold10,"100 Fold"=fold)
    boxplot(dat,main="MSE Distribution",col=c("#4f9da6","#ff5959","#facf5a","#a64f9d","#009E73"),frame=F,cex.axis=1.5,cex.main=1.8,ylim=c(0,0.1))
    text(x=(-0.25),y=0.11, labels=labs[i],cex=4)
}


dev.off()

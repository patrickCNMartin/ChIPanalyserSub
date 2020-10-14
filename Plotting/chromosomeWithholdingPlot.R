pdf("chromosome_withhold_setup.pdf",width=15,height=22)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,8,9,10), ncol=2, byrow=T),height=c(1.4,1.4,1.4,1.4,1.4,1.4,2,2))
par(family="sans",xpd=NA)


for(i in seq_along(trainChIPLoci)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(trainChIPLoci)[i],end(trainChIPLoci)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(trainChIPLoci)[i],(end(trainChIPLoci)[i])+1,by=5000), labels=seq(start(trainChIPLoci)[i],(end(trainChIPLoci)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1) {
  title(main="Training Profiles on chr3R",cex.main=2)
  text(x=start(trainChIPLoci)[i],y=1.05,label="A", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(trainChIPLoci[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(trainChIPLoci[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,trainChIPscores[[i]][seq(0,length(trainChIPscores[[i]]),by=100)],0)
predInd<-c(0,trainPred[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}

for(i in seq_along(validationChIPLociArt)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(validationChIPLociArt)[i],end(validationChIPLociArt)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(validationChIPLociArt)[i],(end(validationChIPLociArt)[i])+1,by=5000), labels=seq(start(validationChIPLociArt)[i],(end(validationChIPLociArt)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1){
  title(main="Validation Profiles on chr3R",cex.main=2)
  text(x=start(validationChIPLociArt)[i],y=1.05,label="B", cex=4)
  text(x=(end(validationChIPLociArt)[i])+2500,y=1.05,label="C", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(validationChIPLociArt[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(validationChIPLociArt[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,validationChIPscoreArt[[i]][seq(0,length(validationChIPscoreArt[[i]]),by=100)],0)
predInd<-c(0,validationPredArt[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}


for(i in seq_along(validationChIPLoci)){
cols<-c("#facf5a","#233142","#ff5959","#4f9da6")


x<-seq(start(validationChIPLoci)[i],end(validationChIPLoci)[i],by=100)
x<-c(x[1]-1,x,x[length(x)]+1)
par(mar=c(6,3,4,1))
plot(0,type="n",axes=F, xlab=' ',ylab=' ',main=' ',xlim=c(head(x,1),tail(x,1)),ylim=c(0,1))
axis(1,at=seq(start(validationChIPLoci)[i],(end(validationChIPLoci)[i])+1,by=5000), labels=seq(start(validationChIPLoci)[i],(end(validationChIPLoci)[i])+1,by=5000),cex.axis=1.6)
#axis(2,at=c(-1.5,-1,-0.5,0),labels=c("Occupancy","Prediction","ChIP","Access"),las=2)
if(i ==1){
  title(main="Validation Profiles on chr2R",cex.main=2)
  text(x=start(validationChIPLoci)[i],y=1.05,label="C", cex=4)
  text(x=(end(validationChIPLoci)[i])+2500,y=1.05,label="C", cex=4)
}
title(xlab=paste0("Genomic Position on : ",seqnames(validationChIPLoci[i])),cex.lab=1.6)


## plotting stuff
noaccess<-.AccessExtract(validationChIPLoci[i],AccessBEAF)[[1]]
for(j in seq_len(nrow(noaccess))){
    rect(noaccess[j,"start"],0,noaccess[j,"end"],0.75,density=50,col=cols[1],lwd=0.8,border=NA)
    #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=135,lwd=0.8,border=NA)
  #rect(noaccess[i,"start"],-0.25,noaccess[i,"end"],0.25,col=cols[1],density=10,angle=90,lwd=0.8,border=NA)
}

chipInd<-c(0,validationChIPscore[[i]][seq(0,length(validationChIPscore[[i]]),by=100)],0)
predInd<-c(0,validationPred[[i]]$ChIP,0)

polygon(x,chipInd,col=cols[2],lwd=2)
lines(x,predInd,col=cols[3],lwd=4)



}

counts<-list(c(1,2),c(3,4),c(5,6),c(7,8))
lims <- list(c(0.5,1),c(-0.5,0.85),c(0.5,1),c(0,0.05))
lablims <-c(1.15,1.05,1.15,0.07)
labs<-LETTERS[4:7]
for(i in seq_along(metricsVal)){
cols <-c("#4f9da6","#95c4c9","#ff5959","#ff9b9b","#facf5a","#fce29c","#009E73","#66c4ab")
#names(metricsVal)<-gsub("\\."," ", names(metricsVal))
#names(metricsVal)<-gsub("MSE","Norm. MSE",names(metricsVal))
par(mar=c(12,8,4,6))
#metricsVal[[i]]<-lapply(metricsVal[[i]])
par(xpd=NA)


boxplot(metricsVal[[i]],main=names(metricsVal)[i],col=cols[counts[[i]]],
frame=F,cex.axis=1.5,cex.main=1.8,ylim=lims[[i]],las=2,cex=2)
title(ylab="Scores", cex.lab=2, line=5)
text(x=0,y=lablims[i], labels=labs[i],cex=4)
#text(x=(-0.55),y=1.1, labels="C",cex=4)

}
dev.off()

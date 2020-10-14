################################################################################
######################           Motif                ##########################
################################################################################



library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)

library(ChIPanalyser)

source("DataHand.R")



## load data

input<-read.table("encode.txt",sep=' ', comment.char='@', stringsAsFactors=F)
pfms<-input[!duplicated(input[,3]),c(3:4)]
pfms<-pfms[c(1:3,16:17),]
name<-strsplit(pfms[,1],"/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/")
for(i in seq_along(name)){
  if(length(name[[i]])==2){
    buffer<-strsplit(name[[i]][2],".pfm")
    print(buffer)
     name[[i]]<-strsplit(name[[i]][2],".pfm")[1]
  }
  else{
    name[[i]]<-name[[i]]
  }
}
#name<-unlist(name)[-c(6,7)]
name<-unlist(name)
#name[c(4,5)]<-c("PC","Pho")
#name<-name[-c(7,8,9,10,11)]
## get the pwm from ChIPanalyser
## drop the unnecssary ones and add csl

#pfms<-pfms[c(1,2,3,4,5,8,14:17),]


name<-c(name,"dfd","CTCF H.Sapiens")
pfms<-rbind(pfms,c("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/pfmDroso/Dfd.pfm","JASPAR"),
                c("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/CTCF.jaspar","JASPAR"))

## dont forget to ad dfd hox if you need to redo this
 pfms <-pfms[!is.na(pfms[,1]),]


pwms<-vector("list",length(name))

for(i in seq_len(nrow(pfms))){
    if(pfms[i,2]=="matrix"){
      pf<-get(load(pfms[i,1]))
    } else{
      pf<-pfms[i,1]
    }

    bufferGPP<-genomicProfiles(PFM=pf,PFMFormat=as.character(pfms[i,2]))
    pwms[[i]]<-PositionFrequencyMatrix(bufferGPP)


}

pwms<-lapply(pwms,function(x){
  (x-min(x))/(max(x)-min(x))
})
#pwms<-pwms[-c(6,7)]
names(pwms)<-name
for(i in seq_along(pwms)){
  local<-apply(pwms[[i]],2,sum)
  if(!all(local==1)){
    diff<-1-local
    for(j in seq_along(diff)){
      if(diff[j]>0){
        pwms[[i]][,j]<-pwms[[i]][,j]+(abs(diff[j])/nrow(pwms[[i]]))
      } else{
        pwms[[i]][,j]<-pwms[[i]][,j]-(abs(diff[j])/nrow(pwms[[i]]))
      }
    }
  }else{
    next
  }
  print(sum(apply(pwms[[i]],2,sum))/ncol(pwms[[i]]))
}







pdf("MotifLogo.pdf")
g<-grid.layout(ncol=3,nrow=2)
for(i in seq_along(pwms)){
  seqLogo(pwms[[i]])
}
dev.off()




mySeqLogo = seqLogo::seqLogo
bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] = NULL
norm = function(x) scale(x, center=FALSE, scale=colSums(x))

pdf("MotifLogo.pdf", width=16, height=10)
grid.newpage()
startxsub5<-0.15
startysub5<-0.47
startxsub10<-0.48
startysub10<-0.47
startxsub15<-0.78
startysub15<-0.47
startxsub20<-0.15
startysub20<- 0.15
for(i in 1:length(pwms)){

    if(i<=2){
      print(i)
      pushViewport(viewport(x=startxsub5, y=startysub5, width=0.34, height=0.33))

      mySeqLogo(pwms[[i]])
      popViewport()
      grid.text(name[i], x=startxsub5+0.03, y=startysub5+0.15, hjust=0.5, vjust=1,gp=gpar(cex=1.5,fontfamilly="mono"))
      startxsub5<-startxsub5
      startysub5<-startysub5+0.35
    }
    if(i>2 & i<=4){
    print(i)
      pushViewport(viewport(x=startxsub10, y=startysub10, width=0.34, height=0.33))

      mySeqLogo(pwms[[i]])
      popViewport()
      grid.text(name[i], x=startxsub10+0.03, y=startysub10+0.15, hjust=0.5, vjust=1,gp=gpar(cex=1.5,fontfamilly="mono"))
      startxsub10<-startxsub10
      startysub10<-startysub10+0.35
    }
    if(i> 4 & i<=6){
      print(i)
      pushViewport(viewport(x=startxsub15, y=startysub15, width=0.34, height=0.33))

      mySeqLogo(pwms[[i]])
      popViewport()
      grid.text(name[i], x=startxsub15+0.03, y=startysub15+0.15, hjust=0.5, vjust=1,gp=gpar(cex=1.5,fontfamilly="mono"))
      startxsub15<-startxsub15
      startysub15<-startysub15+0.35
    }
    if(i==7){
      print(i)
      pushViewport(viewport(x=startxsub20, y=startysub20, width=0.34, height=0.33))

      mySeqLogo(pwms[[i]])
      popViewport()
      grid.text(name[i], x=startxsub20+0.03, y=startysub20+0.15, hjust=0.5, vjust=1,gp=gpar(cex=1.5,fontfamilly="mono"))
      startxsub20<-startxsub20
      startysub20<-startysub20+0.5
    }



}
dev.off()




## A variety of layouts (some a bit mid-bending ...)
layout.torture()
## Demonstration of layout justification
grid.newpage()
testlay <- function(just="centre") {
  pushViewport(viewport(layout=grid.layout(1, 1, widths=unit(1, "inches"),
                          heights=unit(0.25, "npc"),
                          just=just)))
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
  grid.rect()
  grid.text(paste(just, collapse="-"))
  popViewport(2)
}
testlay()
testlay(c("left", "top"))
testlay(c("right", "top"))
testlay(c("right", "bottom"))
testlay(c("left", "bottom"))
testlay(c("left"))
testlay(c("right"))
testlay(c("bottom"))
testlay(c("top"))
# }

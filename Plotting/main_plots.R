###############################################################################
############################# all data parisng ################################
###############################################################################
# the usual although I dont think I actually need this
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


## load that data

aucMatDHS<-get(load("AUCDHS_validated.Rda"))
aucMatcont<-get(load("AUCcont_validated.Rda"))
aucMatNULL<-get(load("AUCNULL_validated.Rda"))


mseMatDHS<-get(load("MSEDHS_validated.Rda"))
mseMatcont<-get(load("MSEcont_validated.Rda"))
mseMatNULL<-get(load("MSENULL_validated.Rda"))

recMatDHS<-get(load("recallDHS_validated.Rda"))
recMatcont<-get(load("recallcont_validated.Rda"))
recMatNULL<-get(load("recallNULL_validated.Rda"))

spearMatDHS<-get(load("spearmanDHS_validated.Rda"))
spearMatcont<-get(load("spearmancont_validated.Rda"))
spearMatNULL<-get(load("spearmanNULL_validated.Rda"))


aucMats<-list("NULL"=aucMatNULL,"Cont"=aucMatcont,"DHS"=aucMatDHS)
mseMats<-list("NULL"=mseMatNULL,"Cont"=mseMatcont,"DHS"=mseMatDHS)
recMats<-list("NULL"=recMatNULL,"Cont"=recMatcont,"DHS"=recMatDHS)
spearMats<-list("NULL"=spearMatNULL,"Cont"=spearMatcont,"DHS"=spearMatDHS)



aucMatDHSTrain<-get(load("AUCDHS_training.Rda"))
aucMatcontTrain<-get(load("AUCcont_training.Rda"))
aucMatNULLTrain<-get(load("AUCNULL_training.Rda"))

mseMatDHSTrain<-get(load("MSEDHS_training.Rda"))
mseMatcontTrain<-get(load("MSEcont_training.Rda"))
mseMatNULLTrain<-get(load("MSENULL_training.Rda"))

recallMatDHSTrain<-get(load("recallDHS_training.Rda"))
recallMatcontTrain<-get(load("recallcont_training.Rda"))
recallMatNULLTrain<-get(load("recallNULL_training.Rda"))

spearMatDHSTrain<-get(load("spearDHS_training.Rda"))
spearMatcontTrain<-get(load("spearcont_training.Rda"))
spearMatNULLTrain<-get(load("spearNULL_training.Rda"))








## functions

ordered<-function(mat,tag=c("CTCF","Hw","BEAF")){
    ord<-c()
    for(i in seq_along(tag)){
       ord<-c(ord,grep(tag[i],colnames(mat)))

    }
    return(mat[,ord])
}


extractMetrics<-function(mat,n=20){

    mat<-mat[seq(10,n+10),]
    return(apply(mat,2, mean))
}

mseOnParam<-function(mat){
   res <-apply(mat,2,function(x){return(x/max(x))})
   res<-apply(mat,2,sd)
   return(res)
}

between<-function(x,low,high,upper){


   if(upper!=high){
      res<- (x>=low & x<high)
   }else{

      res<- (x>low & x<=high)
   }
   return(res)
}

catBuilder<-function(metrics,n=5,dynamic=TRUE,method="AUC"){

   if(dynamic){
      maxi<-max(sapply(metrics,max))
      mini<-min(sapply(metrics,min))
   } else{
      if(method=="AUC"){
      maxi<-1
      mini<-0.35
      }
      if(method=="MSE"){
        maxi<-0.016
        mini<-0.0019
      }
   }
   cats<-seq(mini,maxi, length.out=n+1)

   upper<-cats[n+1]
   buffer<-metrics
   catBuffer<-LETTERS[seq_along(cats)]
   for(i in seq_len(n)){

      buffer[between(buffer,cats[i],cats[i+1],upper)]<-catBuffer[i]
   }
   for(i in seq_along(catBuffer)){
      buffer[buffer==catBuffer[i]]<-i
   }
   name<-names(buffer)
   buffer<-as.numeric(buffer)
   names(buffer)<-name
   if(dynamic){
     return(list(buffer,cats))
   } else{
     return(buffer)
   }
}

formatCoversion<-function(build){
    dhs<-list("DHS"=NULL,"Cont"=NULL,"NULL"=NULL)
    dhs<-lapply(dhs,function(x){rep(0, length(build[[1]][[1]]))})
    for(i in seq_along(build[[1]][[1]])){

       for(j in seq_along(build[[1]])){
          buffer<-c()
          for(k in seq_along(build)){
            buffer<-c(buffer, build[[k]][[j]][i])
          }
          buffer<-median(buffer)
          dhs[[j]][i]<-buffer
       }
    }
   return(dhs)

}

decayExtraction<-function(build,numdata=33){
    bydataset<-vector("list", numdata)
    names(bydataset)<-names(build[[1]][[1]])
    for(i in seq_along(build[[1]][[1]])){
       bydataset[[i]]<-vector("list",3)
       names( bydataset[[i]])<-names(build[[1]])
       for(j in seq_along(build[[1]])){
         for(k in seq_along(build)){
             bydataset[[i]][[j]]<-c(bydataset[[i]][[j]],build[[k]][[j]][[i]])
         }
       }
    }
    return(bydataset)
}

deltaMetric<-function(data,tf=c("CTCF","BEAF","Hw")){
    local<-vector("list", length(data))
    for(i in seq_along(data)){
       local[[i]]<-data[[i]][[3]] - data[[i]][[1]]
    }
    bytf<-vector("list", length(tf))

    for(i in seq_along(tf)){
       bytf[[i]]<-lapply(local,function(x,tf){
                         res<-x[grep(tf,names(x))]
                         return(res)
       },tf[i])
    }

    return(bytf)
}

deltaPval<-function(data){
    pvalMats<-vector("list", length(data))
    ranges<-c(20,50,100,200,500,1000,3273)
    for(i in seq_along(data)){
       pvalMats[[i]]<-matrix(0,ncol=length(data[[i]]),nrow=length(data[[i]]))
       for(j in seq_along(data[[i]])){
         for(k in seq_along(data[[i]])){
            buffer<-t.test(data[[i]][[j]],data[[i]][[k]])$p.value
            if(buffer<0.05){
               pvalMats[[i]][j,k]<-0
            } else{
              pvalMats[[i]][j,k]<-1
            }


         }
       }
       colnames(pvalMats[[i]])<-ranges
       rownames(pvalMats[[i]])<-ranges
    }

     return(pvalMats)
}

overlayMat<-function(mat, method="geometricMean"){
   over<-matrix(0,nrow=20,ncol=17)
   mats<-vector("list",ncol(mat))
   for(i in seq_along(mats)){
     buffer<-mat[,i]
     dim(buffer)<-c(20,17)
     mats[[i]]<-buffer
   }

   for(k in seq_along(mats)){
     buffer<-mats[[k]]
     idx<-seq_along(as.vector(mats[[k]]))


     # define top hits
     top<-head(idx,floor(length(idx)*0.1))

     # Top hits

     if(method %in% c("geometric","MSE","ks","recall")){

         ord<-order(buffer,decreasing=F)
     } else{

         ord<-order(buffer,decreasing=T)
     }
     ##creating a ordered value matrix
     optimalMatrix<-matrix(match(idx,ord),ncol=17, nrow=20)

     ## Location of top hits in ordered value matrix
     topHits<-which(optimalMatrix<max(top),arr.ind=T)

     # Creating empty matrix
     mat<- matrix(0,ncol=17,nrow=20)

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
    rownames(over)<-seq(0.25,5, by =0.25)
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


## getting dat data
#base plots
dhs<-lapply(lapply(aucMats,extractMetrics,n=100),catBuilder,dynamic=F)
#catsDHS<-lapply(dhs,"[[",2)
#dhs<-lapply(dhs,"[[",1)

## with median range
ranges<-c(20,50,100,200,500,1000,3283)

dhsbuild<-vector("list", length(ranges))
for(i in seq_along(ranges)){
    dhsbuild[[i]]<-lapply(lapply(lapply(aucMats,ordered),extractMetrics,n=ranges[i]),catBuilder,dynamic=F,method="AUC")

}
tags<-names(dhsbuild[[1]][[1]])
dhs<-formatCoversion(dhsbuild)
dhs<-lapply(dhs,function(x,tag){names(x)<-tag;return(x)},tags)

### Loading data from hat maps

mse<-lapply(lapply(lapply(mseMats,ordered),mseOnParam),catBuilder,dynamic=F,method="MSE")








# for heat maps


### plot stuff i guess

tfs<-c("CTCF","Hw","BEAF")

n=5
ptCol<-seq_len(n)
ptSize<-seq(4,12,length.out=n)

colfunc<-colorRampPalette(c("#facf5a","#ff5959","#233142"))
cols<-colfunc(5)

tags<-names(dhsbuild[[1]][[1]])
tags<-gsub("_"," ",tags)
tags<-gsub("modEncode","",tags,ignore.case=T)
tags<-gsub("Su","su",tags)
catsDHS<-seq(0.43,1, length.out=6)
catsMSE<-seq(0.0019,0.016,length.out=6)
labs<-LETTERS[1:3]

pdf("AccessMapsAllRegions.pdf",width=19,height=13.5)
layout(matrix(c(rep(1,6),rep(2,6),5,rep(3,5),6,4,4,4,4,4),byrow=F,ncol=4),
width=c(1,1,1,1))

par(family="sans")
par(xpd=NA)



for(i in seq_along(tfs)){
    tags<-names(dhs[[1]])
    localTF<-lapply(dhs,function(dhs,tf){dhs[grep(tf,names(dhs))]},tfs[i])
    localMSE<-lapply(mse,function(mse,tf){mse[grep(tf,names(mse))]},tfs[i])
    localTF<-do.call("cbind", localTF)
    localMSE<-do.call("cbind",localMSE)
    tags<-rownames(localTF)
    tags<-gsub("_"," ",tags)
    tags<-gsub("modEncode","",tags,ignore.case=T)
    tags<-gsub("Su","su",tags)
    if(i !=3){
      par(mar=c(20,18.5,8,4))
    }else{
      par(mar=c(20,18.5,0,4))
    }#

    plot(0,type="n",ylim=c(0.5,nrow(localTF)+1),xlim=c(0.4,3),axes=F,xlab="",ylab="")
    axis(2,at=seq_len(nrow(localTF)),labels=tags,las=2,cex.axis=1.6)
    axis(1,at=c(1,2,3),labels=F,line=(-1))
    y<-(-0.02*nrow(localTF))
    text(x=c(1,2,3),y =y, srt = 45, adj = 1,labels = rev(c("DHS Accessibility","Continuous Accessibility","No Accessibility")),
    xpd = TRUE,cex=2.2)
    text(x=-1, y=nrow(localTF)+0.1*(nrow(localTF)), labels=labs[i],cex=5,xpd=T)

    for(j in seq_len(nrow(localTF))){

        points(x=1,y=j,cex=ptSize[localMSE[j,1]],col=cols[localTF[j,1]],pch=19)
        points(x=2,y=j,cex=ptSize[localMSE[j,2]],col=cols[localTF[j,2]],pch=19)
        points(x=3,y=j,cex=ptSize[localMSE[j,3]],col=cols[localTF[j,3]],pch=19)
    }

}
par(mar=c(14,16,32,14))
par(xpd=NA)
plot(0,type="n",ylim=c(0,-6),xlim=c(0,1),axes=F,xlab="",ylab="")
axis(4,at=seq(-1,-6),labels=round(catsDHS,digits=3),las=2,cex.axis=2)
rect(0,-6,1,-5,col=cols[5])
rect(0,-5,1,-4,col=cols[4])
rect(0,-4,1,-3,col=cols[3])
rect(0,-3,1,-2,col=cols[2])
rect(0,-2,1,-1,col=cols[1])
text(x=0.5,y=-0.5,"Median AUC over Validation",cex=2)

par(mar=c(4,4,4,4))
plot(0,type="n",axes=F,xlab="",ylab="")

par(mar=c(9.5,6,2,2))
par(xpd=NA)
plot(0,type="n",ylim=c(0,1),xlim=c(0,1.2),axes=F,xlab="",ylab="")
axis(4,at=seq(-5,-1,by=1),labels=paste(round(catsMSE[1:5],digits=4),"-",round(catsMSE[2:6],digits=4)),las=2,cex.axis=1.8,tick=F,line=(-11))
points(x=0.55,y=-5,cex=ptSize[1])
points(x=0.55,y=-4,cex=ptSize[2])
points(x=0.55,y=-3,cex=ptSize[3])
points(x=0.55,y=-2,cex=ptSize[4])
points(x=0.55,y=-1,cex=ptSize[5])
text(y=0.1,x=0.55,"SD MSE over Training",cex=2)
dev.off()


################################################################################
############################ Decay box plots ###################################
################################################################################

ranges<-c(20,50,100,200,500,1000,3283)

dhsbuild<-vector("list", length(ranges))
for(i in seq_along(ranges)){
    dhsbuild[[i]]<-lapply(lapply(aucMats,ordered),extractMetrics,n=ranges[i])

}

dhs<-deltaMetric(dhsbuild)


dhsPval<-deltaPval(dhs)



pdf("boxplotDecayAUCpval.pdf",width=18,height=16)
par(family="sans")
par(xpd=NA)
layout(matrix(1:6,ncol=2,byrow=F),width=c(4,1.6))
par(mar=c(8,8,8,4))
mains<-c("CTCF","BEAF-32","su(Hw)")
cols<-c("#facf5a","#ff5959","#4f9da6")
tagsauc<-LETTERS[1:3]
for(i in seq_along(dhs)){
  names(dhs[[i]])<-c(20,50,100,200,500,1000,3283)
  boxplot(dhs[[i]],las=2,col=cols[i],frame=F,cex.axis=1.7,cex.names=1.4,ylim=c(-0.5,0.3))
  lines(y=c(0,0),x=c(0.3,7.4),lty=2,lwd=1.8)
  title(xlab="Number of Regions Selected", line=6,cex.lab=1.8)
  title(ylab="Delta mean AUC Score (DHS - NULL)", line=6,cex.lab=1.8)
  title(main=mains[i],cex.main=1.6)

  text(x=0,y=0.5,tagsauc[i],cex=4)

}
tagsauc<-LETTERS[4:6]
colfunc<-colorRampPalette(c("#4f9da6","#edf5f6"))
cols<-colfunc(30)
par(mar=c(8,2,4,4))
for(i in seq_along(dhsPval)){

    image(x=1:7,y=1:7,t(dhsPval[[i]]),zlim=c(0,1),col=cols,axes=F,xlab="",ylab="")
    axis(1,at=1:7,labels=c(20,50,100,200,500,1000,3283),cex.axis=1.8)
    axis(2,at=1:7,labels=c(20,50,100,200,500,1000,3283),cex.axis=1.8,las=2)
    text(x=-0.8,y=7.8,tagsauc[i],cex=4)
    title(main="Differences between # regions - (p<0.05) ",cex.main=1.8)

}

dev.off()


################################################################################
############################ Data consistency ##################################
################################################################################

trainMats<-list("NULL"=aucMatNULLTrain,"Cont"=aucMatcontTrain,"DHS"=aucMatDHSTrain)
trainMats<-list("NULL"=mseMatNULLTrain,"Cont"=mseMatcontTrain,"DHS"=mseMatDHSTrain)
trainMats<-list("NULL"=recallMatNULLTrain,"Cont"=recallMatcontTrain,"DHS"=recallMatDHSTrain)
trainMats<-list("NULL"=spearMatNULLTrain,"Cont"=spearMatcontTrain,"DHS"=spearMatDHSTrain)


overlay<-lapply(trainMats,bytfmat)

mains<-lapply(overlay, names)


## dont forget to change in the lappy thingy the method you are using for extrcation
overlay<-lapply(1:3,function(x,mats){
                          local<-mats[[x]]
                          idx<-seq_along(local)
                          over<-lapply(idx,function(idx,local){
                                     return(overlayMat(local[[idx]],method="spearman"))},local)

                          return(over)},overlay)



colfunc<-colorRampPalette(c("white","#ff5959","#233142"))
labs<-LETTERS[1:9]

pdf("overlayMatriciesSpearman.pdf",height=18,width=18)
par(family="sans")
par(xpd=NA)
layout(matrix(1:18,ncol=6,byrow=T),width=rep(c(1,0.28),times=3))
xlabs<-c(1, 10, 20, 50, 100,
    200, 500,1000,2000, 5000,10000,20000,50000, 100000,
    200000, 500000, 1000000)
ylabs<-seq(0.25,5,by=0.25)
for(j in seq_along(overlay)){
  for(i in seq_along(overlay[[j]])){
    par(mar=c(8,10,8,0))
  Colors <- colfunc(15)
  legend_image<-as.raster(matrix(rev(Colors),ncol=1))
  graphics::image(1:length(xlabs),1:length(ylabs),t(overlay[[j]][[i]]),
  axes = FALSE, xlab=" ", ylab=" ",col=Colors)
  text(x=(-1), y=length(ylabs)+2,labels=labs[i],cex=5)
  mainTit<-gsub("Hw","su(Hw)",mains[[j]])
  title(main=mainTit[i],cex.main=2)
  title( ylab="Scaling Factor", line=6, cex.lab=1.6)
  title(xlab="Number of Bound Molecules",line=6.5, cex.lab=1.6)
  axis(1,at=seq_along(xlabs),labels=F)
  text(seq_along(xlabs),y =(-0.2), srt = 45, adj = 1,labels = xlabs, xpd = TRUE,cex=1.4)
  axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,las= HORIZONTAL<-1,cex.axis=1.4)
  par(mar=c(5,1,6,0))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  text(x=1.6, y =seq(0,1,l=5) , labels = round(seq(min(overlay[[j]][[i]]),max(overlay[[j]][[i]]),l=5),2),cex=1.4)
  rasterImage(legend_image, 0, 0, 1,1)
  }
}
dev.off()


### decay curves for all of them

aucMats<-list("NULL"=aucMatNULL,"Cont"=aucMatcont,"DHS"=aucMatDHS)
mseMats<-list("NULL"=mseMatNULL,"Cont"=mseMatcont,"DHS"=mseMatDHS)
recMats<-list("NULL"=recMatNULL,"Cont"=recMatcont,"DHS"=recMatDHS)
spearMats<-list("NULL"=spearMatNULL,"Cont"=spearMatcont,"DHS"=spearMatDHS)

ranges<-c(20,50,100,200,500,1000,3283)

mats<-recMats

dhsbuild<-vector("list", length(ranges))
for(i in seq_along(ranges)){
    dhsbuild[[i]]<-lapply(lapply(mats,ordered),extractMetrics,n=ranges[i])

}
decay<-decayExtraction(dhsbuild)

pdf("decaySpearman.pdf",width=32, height=34)
cols<-c("#233142","#ff5959","#4f9da6")
par(family="sans")
par(xpd=NA)
mat<-matrix(c(1:14,0,15:25,0,0,0,0,26:33,0,0),ncol=5,nrow=8, byrow=T)
mat[which(mat==0,arr.ind=T)]<-34:40
layout(mat)
ranges<-c(20,50,100,200,500,1000,3283)
count<-1
labels<-LETTERS[1:3]
for(i in seq_along(decay)){
   par(mar=c(5.5,9,10,11))
   ## for geo and mse
   #ylims<-c(min(min(log(decay[[i]]$medDHS)),min(log(decay[[i]]$medNULL))),
          #max(max(log(decay[[i]]$medDHS)),max(log(decay[[i]]$medNULL))))
  # for AUC
  #ylims<-c(min(min(decay[[i]]$medDHS),min(decay[[i]]$medNULL)),
         #max(max(decay[[i]]$medDHS),max(decay[[i]]$medNULL)))
   #ylims<-log(c(0.0008,0.06)+1)
   ylims<-c(0.3,1)
   xlims<-c(1,length(ranges))

   plot(0,type="n",axes=F,xlab="",ylab="",xlim=xlims,
   ylim=ylims)
   if(i %in% c(1,15,26)){text(x=-0.6,y=ylims[2]+0.25*ylims[2],labels[count],cex=7);count<-count+1}
   #if(i %in% c(1,15,26)){text(x=-0.6,y=-2,labels[count],cex=7);count<-count+1}

   main<-names(decay)[[i]]
   main<-gsub("_"," ",main)
   main<-gsub("modEncode ","",main,ignore.case=T)

   main<-gsub("Su","su",main)
   axis(1,at=seq_len(xlims[2]),labels=ranges,cex.axis=1.6)
   axis(2,at=round(seq(ylims[1],ylims[2],length.out=5),digits=2),labels=round(seq(ylims[1],ylims[2],length.out=5),digits=2),
   cex.axis=1.6,las=2)
   title(xlab="# of Loci",cex.lab=1.8)
   title(ylab="max Spearman score",cex.lab=1.8,line=6)
   title(main=main,cex.main=2)
   ## for geo and mse
   #lines(seq_along(decay[[i]]$medDHS),log(decay[[i]]$medDHS),col=cols[1],type="b",lwd=2)
   #lines(seq_along(decay[[i]]$medNULL),log(decay[[i]]$medNULL),col=cols[2],type="b",lwd=2)
   ##for auc
   lines(seq_along(decay[[i]][[1]]),decay[[i]][[1]],col=cols[1],type="b",lwd=2)
   lines(seq_along(decay[[i]][[2]]),decay[[i]][[2]],col=cols[2],type="b",lwd=2)
   lines(seq_along(decay[[i]][[3]]),decay[[i]][[3]],col=cols[3],type="b",lwd=2)

    legY<-ylims[1]+((ylims[2]-ylims[1])/2)
    print(legY)
   legend(x=xlims[2],y=legY,legend=c("No Access","Cont. Access","DHS Access"),fill=c(cols[1],cols[2],cols[3]),bty="n",cex=1.6)
}

par(xpd=NA)
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),xlim=c(0,1),ylim=c(1,2),type="l",lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),xlim=c(0,1),ylim=c(1,2),type="l",lty=2,lwd=2,axes=F,xlab="",ylab="")
dev.off()


#### let's do some data adding shall
